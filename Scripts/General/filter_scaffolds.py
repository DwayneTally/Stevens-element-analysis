#!/usr/bin/env python3
"""
For each flagged genome assembly, keeps only the top n scaffolds by length
(chromosome number n), discarding unplaced scaffolds and Y chromosomes.

Output files are written alongside the originals with a .filtered.fna suffix,
unless an output directory is specified with --outdir.

Usage:
    python filter_scaffolds.py --indir /path/to/fastas --outdir /path/to/output
    python filter_scaffolds.py --indir /path/to/fastas --dry-run
"""

import os
import sys
import argparse

# ---------------------------------------------------------------------------
# Assemblies to filter: filename_stem -> n chromosomes to keep
# Bruchidius siliquastri is listed but already has a _clean version —
# the script will skip it if a .filtered or _clean file already exists.
# ---------------------------------------------------------------------------
TO_FILTER = {
    # Chrysomeloidea — large excesses
    "Altica_lythri_GCA_964028335.1_genomic":              12,
    "Arhopalus_rusticus_GCA_964213965.1_genomic":          9,
    "Chrysomela_aeneicollis_GCA_027562985.1_genomic":     21,
    "Crioceris_asparagi_GCA_958507055.1_genomic":          8,
    "Callosobruchus_maculatus_GCA_040182625.1_genomic":    9,
    "Octodonta_nipae_GCA_034190945.1_genomic":             9,
    "Agelastica_alni_GCA_950111635.2_genomic":            12,
    "Lochmaea_crataegi_GCA_947563755.1_genomic":          16,
    "Cryptocephalus_primarius_GCA_963576515.1_genomic":   20,
    "Leptinotarsa_decemlineata_GCA_024712935.1_genomic":  18,
    "Rutpela_maculata_GCA_936432065.2_genomic":           10,
    "Stenurella_melanura_GCA_963583905.1_genomic":        10,
    "Stictoleptura_scutellata_GCA_964212005.1_genomic":   10,
    "Leptura_quadrifasciata_GCA_963675555.1_genomic":     10,
    "Rhagium_mordax_GCA_963680705.1_genomic":             10,
    "Diorhabda_carinulata_GCF_026250575.1_genomic":       14,
    "Chrysolina_haemoptera_GCA_958298965.1_genomic":      20,
    "Crepidodera_aurea_GCA_949320105.2_genomic":          10,
    "Phyllotreta_striolata_GCA_918026865.1_genomic":      15,
    "Brontispa_longissima_GCA_040580785.1_genomic":       10,
    "Galeruca_laticollis_GCA_963921935.1_genomic":        12,
    "Ophraella_communa_GCA_035357415.1_genomic":          17,
    "Galerucella_nymphaeae_GCA_963978555.1_genomic":      16,
    "Plagiodera_versicolora_GCA_035321685.1_genomic":     16,
    "Pogonocherus_hispidulus_GCA_963924545.1_genomic":    10,
    "Psylliodes_chrysocephala_GCA_927349885.1_genomic":   23,  # note: filename uses chrysocephala
    "Gastrophysa_polygoni_GCA_963576655.1_genomic":       12,
    "Monochamus_alternatus_GCA_037114965.1_genomic":      10,
    "Tetropium_fuscum_GCA_964058775.1_genomic":           12,
    "Phaedon_cochleariae_GCA_918026855.4_genomic":        17,
    "Psylliodes_chrysocephalus_GCA_927349885.1_genomic":  23,  # fallback spelling
    "Bruchidius_siliquastri_GCA_949316355.1_genomic":     10,
}


def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file. Memory-efficient."""
    with open(path, "r") as fh:
        header = None
        parts = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line
                parts = []
            else:
                parts.append(line)
        if header is not None:
            yield header, "".join(parts)


def write_fasta(records, path, line_width=60):
    """Write (header, seq) records to a FASTA file with wrapped sequences."""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")


def filter_fasta(in_path, out_path, keep_n, dry_run=False):
    print(f"  Reading  {os.path.basename(in_path)} ...", end=" ", flush=True)
    records = list(parse_fasta(in_path))
    records.sort(key=lambda r: len(r[1]), reverse=True)
    total = len(records)
    kept = records[:keep_n]
    dropped = total - keep_n
    print(f"{total} scaffolds → keeping {keep_n}, dropping {dropped}")
    if dry_run:
        for i, (h, s) in enumerate(kept):
            name = h.split()[0][1:]
            print(f"    [{i+1:>3}] {name}  ({len(s):,} bp)")
        return
    write_fasta(kept, out_path)
    print(f"  Written  {os.path.basename(out_path)}")


def main():
    parser = argparse.ArgumentParser(description="Filter genome FASTAs to top-n scaffolds by length.")
    parser.add_argument("--indir", required=True, help="Directory containing input .fna files")
    parser.add_argument("--outdir", default=None, help="Output directory (default: same as indir)")
    parser.add_argument("--dry-run", action="store_true", help="Print what would be done without writing files")
    parser.add_argument("--suffix", default=".filtered.fna", help="Output file suffix (default: .filtered.fna)")
    args = parser.parse_args()

    indir = args.indir.rstrip("/")
    outdir = args.outdir.rstrip("/") if args.outdir else indir
    os.makedirs(outdir, exist_ok=True)

    skipped = []
    processed = []
    not_found = []

    for stem, n in sorted(TO_FILTER.items()):
        in_path = os.path.join(indir, stem + ".fna")
        out_path = os.path.join(outdir, stem + args.suffix)

        # Skip if a filtered/clean version already exists
        already_done = [
            os.path.join(outdir, stem + ".filtered.fna"),
            os.path.join(indir, stem + "_clean.fna"),
            os.path.join(outdir, stem + "_clean.fna"),
        ]
        if any(os.path.exists(p) for p in already_done) and not args.dry_run:
            print(f"  Skipping {stem} — filtered file already exists")
            skipped.append(stem)
            continue

        if not os.path.exists(in_path):
            print(f"  MISSING  {in_path}")
            not_found.append(stem)
            continue

        print(f"\n{stem}  (keep top {n})")
        filter_fasta(in_path, out_path, n, dry_run=args.dry_run)
        processed.append(stem)

    print(f"\n{'='*60}")
    print(f"Done.  Processed: {len(processed)}  Skipped (already done): {len(skipped)}  Not found: {len(not_found)}")
    if not_found:
        print("Missing files:")
        for s in not_found:
            print(f"  {s}.fna")


if __name__ == "__main__":
    main()
