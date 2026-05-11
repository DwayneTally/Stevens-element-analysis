#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def norm_name(name):
    name = name.strip()
    if name.endswith(".bed"):
        name = name[:-4]
    if name.endswith("_genomic"):
        name = name[:-8]
    return name


def read_exclude_list(path):
    if path is None:
        return set()
    with open(path, "r", encoding="utf-8") as f:
        return {norm_name(line) for line in f if line.strip()}


def build_earlgrey_map(base_dir):
    mapping = {}
    for genome_dir in base_dir.iterdir():
        if not genome_dir.is_dir():
            continue
        summary_dir = genome_dir / "summaryFiles"
        if not summary_dir.exists():
            continue
        matches = list(summary_dir.glob("*.softmasked.fasta"))
        if matches:
            mapping[norm_name(genome_dir.name)] = matches[0]
    return mapping


def build_bed_map(bed_dir):
    mapping = {}
    for bed in bed_dir.glob("*.bed"):
        mapping[norm_name(bed.name)] = bed
    return mapping


def read_x_scaffolds_from_bed(bed_file):
    scaffolds = set()
    with open(bed_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts:
                scaffolds.add(parts[0])
    return scaffolds


def get_fasta_headers(fasta_file):
    headers = []
    with open(fasta_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def get_species_specific_x_scaffolds(genome, fasta_headers, bed_x_scaffolds):
    """
    Return the set of scaffolds to treat as X for this genome.
    Default: use the BED scaffolds.

    Special cases:
    1. Crioceris_asparagi_GCA_958507055.1
       Only consider chromosome-level OY scaffolds at all.
       X remains the BED scaffold(s), but autosomes will be limited to OY too.

    2. Cryptophagus_acutangulus_GCA_963556235.1
       BED X scaffolds are missing from FASTA.
       Use chromosome-like OY scaffolds and infer X as whatever is left
       after treating the first 13 OY scaffolds as autosomes.
       This assumes 13 autosomes and 2 X chromosomes.
    """
    fasta_header_set = set(fasta_headers)

    if genome == "Cryptophagus_acutangulus_GCA_963556235.1":
        oy_headers = sorted([h for h in fasta_headers if h.startswith("OY")])

        # If BED X scaffolds are missing, infer X as the remaining OY scaffolds
        # after the autosomes. Based on your observation: correct autosome number,
        # remaining chromosome-level scaffolds should be X.
        if not (bed_x_scaffolds & fasta_header_set):
            if len(oy_headers) >= 15:
                autosomes = set(oy_headers[:13])
                inferred_x = set(oy_headers) - autosomes
                print(f"INFO: {genome} BED X scaffolds not found in FASTA; inferring X as remaining OY scaffolds: {sorted(inferred_x)}")
                return inferred_x

    return bed_x_scaffolds


def fasta_header_allowed(genome, header):
    """
    Species-specific filtering of which FASTA headers to include at all.
    """
    if genome == "Crioceris_asparagi_GCA_958507055.1":
        return header.startswith("OY")
    if genome == "Cryptophagus_acutangulus_GCA_963556235.1":
        return header.startswith("OY")
    return True


def count_masking_by_compartment(fasta_file, x_scaffolds, genome):
    counts = {
        "X_noNeoX": {"masked": 0, "unmasked": 0, "seqs": 0},
        "Autosomes": {"masked": 0, "unmasked": 0, "seqs": 0},
    }

    current_group = None
    current_allowed = False

    with open(fasta_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                current_allowed = fasta_header_allowed(genome, header)

                if not current_allowed:
                    current_group = None
                    continue

                current_group = "X_noNeoX" if header in x_scaffolds else "Autosomes"
                counts[current_group]["seqs"] += 1
                continue

            if current_group is None:
                continue

            seq = line.strip()
            for ch in seq:
                if ch.isalpha():
                    if ch.islower():
                        counts[current_group]["masked"] += 1
                    elif ch.isupper():
                        counts[current_group]["unmasked"] += 1

    return counts


def calc_masked_percent(masked, unmasked):
    total = masked + unmasked
    if total == 0:
        return 0, 0.0
    return total, round((masked / total) * 100, 2)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate soft-masked percent for ancestral X (no neoX) and autosomes using BED files."
    )
    parser.add_argument(
        "--earlgrey-base",
        required=True,
        help="Base directory containing EarlGrey genome folders with summaryFiles/*.softmasked.fasta"
    )
    parser.add_argument(
        "--x-bed-dir",
        required=True,
        help="Directory containing BED files for ancestral X chromosomes"
    )
    parser.add_argument(
        "--exclude-list",
        default=None,
        help="Optional text file of genomes to exclude (e.g. neosex.txt)"
    )
    parser.add_argument(
        "--out",
        default="x_vs_autosomes_softmasked.csv",
        help="Output CSV file"
    )
    args = parser.parse_args()

    earlgrey_base = Path(args.earlgrey_base)
    x_bed_dir = Path(args.x_bed_dir)

    if not earlgrey_base.exists():
        raise FileNotFoundError(f"EarlGrey base directory not found: {earlgrey_base}")
    if not x_bed_dir.exists():
        raise FileNotFoundError(f"X BED directory not found: {x_bed_dir}")

    earlgrey_map = build_earlgrey_map(earlgrey_base)
    bed_map = build_bed_map(x_bed_dir)
    exclude = read_exclude_list(args.exclude_list)

    output_rows = []
    missing_softmasked = []
    excluded = []

    genomes_to_process = sorted(bed_map.keys())

    for genome in genomes_to_process:
        if genome in exclude:
            excluded.append(genome)
            continue

        fasta = earlgrey_map.get(genome)
        bed = bed_map.get(genome)

        if fasta is None:
            missing_softmasked.append(genome)
            continue

        bed_x_scaffolds = read_x_scaffolds_from_bed(bed)
        if not bed_x_scaffolds:
            print(f"WARNING: no scaffolds found in BED for {genome}")
            continue

        fasta_headers = get_fasta_headers(fasta)
        x_scaffolds = get_species_specific_x_scaffolds(genome, fasta_headers, bed_x_scaffolds)

        matched_x = sorted(set(x_scaffolds) & set(fasta_headers))
        if len(matched_x) == 0:
            print(f"WARNING: no X scaffolds matched FASTA for {genome}")
            print(f"  BED/inferred X scaffolds: {sorted(x_scaffolds)}")
            continue

        counts = count_masking_by_compartment(fasta, x_scaffolds, genome)

        for group in ["X_noNeoX", "Autosomes"]:
            masked = counts[group]["masked"]
            unmasked = counts[group]["unmasked"]
            total_bp, masked_percent = calc_masked_percent(masked, unmasked)

            output_rows.append({
                "genome": genome,
                "group": group,
                "masked_bp": masked,
                "unmasked_bp": unmasked,
                "total_bp": total_bp,
                "masked_percent": masked_percent,
                "n_fasta_records": counts[group]["seqs"],
                "matched_x_scaffolds": ",".join(matched_x),
                "softmasked_fasta": str(fasta),
                "x_bed": str(bed),
            })

    with open(args.out, "w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=[
                "genome",
                "group",
                "masked_bp",
                "unmasked_bp",
                "total_bp",
                "masked_percent",
                "n_fasta_records",
                "matched_x_scaffolds",
                "softmasked_fasta",
                "x_bed",
            ],
        )
        writer.writeheader()
        writer.writerows(output_rows)

    processed = sorted({r["genome"] for r in output_rows})

    print(f"Wrote {len(output_rows)} rows for {len(processed)} genomes to {args.out}")

    if processed:
        print("\nProcessed genomes:")
        for g in processed:
            print(f"  {g}")

    if missing_softmasked:
        print("\nMissing softmasked FASTA:")
        for g in missing_softmasked:
            print(f"  {g}")

    if excluded:
        print("\nExcluded by exclude-list:")
        for g in excluded:
            print(f"  {g}")


if __name__ == "__main__":
    main()
