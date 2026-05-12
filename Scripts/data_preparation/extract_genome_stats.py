#!/usr/bin/env python3
import sys
import re
import csv
import gzip
from pathlib import Path

BASE_DIR = Path(".").resolve()
FASTA_DIR = (BASE_DIR / "../unique_fna_clean").resolve()

GFF_DIR_SUFFIX = "_GFF3_filtered_dropIDs"
FA_EXTS = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")

ACC_RE = re.compile(r"(GCA|GCF)_\d+(?:\.\d+)?")

def open_text_maybe_gz(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")

def fasta_stats(fa_path: Path):
    """
    chromosomes = number of FASTA records
    genome_bp   = sum of sequence lengths
    gc_percent  = 100 * (G+C) / (A+C+G+T)   (ignores N/other)
    """
    chrom = 0
    genome_bp = 0
    gc = 0
    acgt = 0

    with open_text_maybe_gz(fa_path) as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                chrom += 1
                continue
            seq = line.strip().upper()
            if not seq:
                continue
            genome_bp += len(seq)
            # fast counting without per-char loops
            g = seq.count("G")
            c = seq.count("C")
            a = seq.count("A")
            t = seq.count("T")
            gc += (g + c)
            acgt += (a + c + g + t)

    gc_percent = (100.0 * gc / acgt) if acgt else 0.0
    return chrom, genome_bp, gc_percent

def gff_gene_count(gff_path: Path):
    """
    Prefer 'gene' features. If none, fall back to 'mRNA'.
    """
    gene = 0
    mrna = 0
    with gff_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            t = parts[2]
            if t == "gene":
                gene += 1
            elif t == "mRNA":
                mrna += 1
    return gene if gene > 0 else mrna

def find_fasta_for_gff(gff3_path: Path):
    """
    Match FASTA primarily by GCA/GCF accession in filename.
    If not present, try a genus_species prefix fallback.
    """
    gname = gff3_path.name

    # 1) Accession match (best)
    m = ACC_RE.search(gname)
    if m:
        acc = m.group(0)
        hits = []
        for p in FASTA_DIR.iterdir():
            if not p.is_file():
                continue
            lname = p.name.lower()
            if acc.lower() in lname and any(lname.endswith(ext) for ext in FA_EXTS):
                hits.append(p)
        if hits:
            # prefer unzipped if both exist
            hits_sorted = sorted(hits, key=lambda x: (x.name.lower().endswith(".gz"), len(x.name)))
            return hits_sorted[0]

    # 2) Fallback: genus_species prefix
    # e.g., "Tribolium_castaneum_..." -> take first two underscore tokens
    toks = gname.split("_")
    if len(toks) >= 2:
        prefix = f"{toks[0]}_{toks[1]}".lower()
        hits = []
        for p in FASTA_DIR.iterdir():
            if not p.is_file():
                continue
            lname = p.name.lower()
            if lname.startswith(prefix) and any(lname.endswith(ext) for ext in FA_EXTS):
                hits.append(p)
        if hits:
            hits_sorted = sorted(hits, key=lambda x: (x.name.lower().endswith(".gz"), len(x.name)))
            return hits_sorted[0]

    return None

def iter_gff3_files(dir_path: Path):
    for p in sorted(dir_path.iterdir()):
        if p.is_file() and p.name.lower().endswith((".gff3", ".gff")):
            yield p

def main():
    if not FASTA_DIR.exists():
        print(f"[ERROR] FASTA_DIR not found: {FASTA_DIR}", file=sys.stderr)
        sys.exit(1)

    out_csv = Path("genome_stats_by_species.csv")
    rows = []
    missing_fasta = 0
    processed = 0

    for d in sorted(BASE_DIR.iterdir()):
        if not d.is_dir() or not d.name.endswith(GFF_DIR_SUFFIX):
            continue

        group = d.name[: -len(GFF_DIR_SUFFIX)]  # superfamily/family label from directory name

        for gff3 in iter_gff3_files(d):
            processed += 1

            fa = find_fasta_for_gff(gff3)
            if fa is None:
                missing_fasta += 1
                print(f"[WARN] No FASTA match for GFF3: {gff3} (group={group})", file=sys.stderr)
                continue

            try:
                chrom, genome_bp, gc_percent = fasta_stats(fa)
            except Exception as e:
                print(f"[WARN] Failed FASTA parse: {fa} ({e})", file=sys.stderr)
                continue

            genes = gff_gene_count(gff3)
            genome_mb = genome_bp / 1_000_000.0

            rows.append({
                "group": group,
                "species_file": gff3.name,
                "chromosomes": chrom,
                "genome_size_mb": genome_mb,
                "gc_percent": gc_percent,
                "genes": genes,
                "fasta": str(fa),
                "gff3": str(gff3),
            })

    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "group","species_file","chromosomes","genome_size_mb","gc_percent","genes","fasta","gff3"
        ])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {out_csv} with {len(rows)} rows", file=sys.stderr)
    print(f"Processed {processed} GFF3 files; missing FASTA matches for {missing_fasta}", file=sys.stderr)

if __name__ == "__main__":
    main()

