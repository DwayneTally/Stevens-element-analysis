#!/usr/bin/env python3

import os
import shutil
import argparse
import pandas as pd


def parse_fasta(path):
    with open(path, "r") as fh:
        header = None
        seq = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)


def write_fasta(records, outpath, width=60):
    with open(outpath, "w") as out:
        for header, seq in records:
            out.write(header + "\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")


def filter_fasta(infile, outfile, keep_n, dry_run=False):
    records = list(parse_fasta(infile))
    records.sort(key=lambda x: len(x[1]), reverse=True)

    total = len(records)
    kept = records[:keep_n]

    total_bp = sum(len(r[1]) for r in records)
    kept_bp = sum(len(r[1]) for r in kept)

    print(f"FILTER: {os.path.basename(infile)}")
    print(f"  scaffolds: {total} -> {keep_n}")
    print(f"  bp: {total_bp:,} -> {kept_bp:,}")

    if not dry_run:
        write_fasta(kept, outfile)


def safe_copy(src, dst, dry_run=False):
    print(f"COPY: {os.path.basename(src)} -> {os.path.basename(dst)}")
    if not dry_run:
        shutil.copy2(src, dst)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True)
    parser.add_argument("--indir", required=True,
                        help="Directory containing original and/or cleaned FASTAs")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--missing-log", default="missing_genomes.txt")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_excel(args.table)

    species_col = "Species"
    chrom_col = "Chromosome number (n)"

    all_fastas = [
        f for f in os.listdir(args.indir)
        if f.endswith(".fna")
    ]

    processed = 0
    copied = 0
    skipped_existing = 0
    missing = 0
    no_chrom = 0

    missing_lines = []

    for _, row in df.iterrows():
        species = str(row[species_col]).strip()

        if pd.isna(row[chrom_col]):
            print(f"SKIP no chromosome count: {species}")
            no_chrom += 1
            continue

        keep_n = int(row[chrom_col])
        species_key = species.replace(" ", "_")

        matches = [f for f in all_fastas if species_key in f]

        if not matches:
            print(f"MISSING: {species}")
            missing_lines.append(species)
            missing += 1
            continue

        # Prefer files that are already cleaned/filtered
        clean_matches = [
            f for f in matches
            if "_clean" in f or "filtered" in f
        ]

        raw_matches = [
            f for f in matches
            if "_clean" not in f and "filtered" not in f
        ]

        if clean_matches:
            infile_name = clean_matches[0]
            infile = os.path.join(args.indir, infile_name)
            outfile = os.path.join(args.outdir, infile_name)

            if os.path.exists(outfile):
                print(f"SKIP already in outdir: {infile_name}")
                skipped_existing += 1
                continue

            safe_copy(infile, outfile, args.dry_run)
            copied += 1
            continue

        if raw_matches:
            infile_name = raw_matches[0]
            infile = os.path.join(args.indir, infile_name)
            outfile_name = infile_name.replace(".fna", ".filtered.fna")
            outfile = os.path.join(args.outdir, outfile_name)

            if os.path.exists(outfile):
                print(f"SKIP already in outdir: {outfile_name}")
                skipped_existing += 1
                continue

            filter_fasta(infile, outfile, keep_n, args.dry_run)
            processed += 1
            continue

    with open(args.missing_log, "w") as out:
        for species in missing_lines:
            out.write(species + "\n")

    print("\nSUMMARY")
    print(f"Filtered new genomes: {processed}")
    print(f"Copied clean/filtered genomes: {copied}")
    print(f"Skipped already in outdir: {skipped_existing}")
    print(f"Missing FASTA: {missing}")
    print(f"No chromosome count: {no_chrom}")
    print(f"Missing log written to: {args.missing_log}")


if __name__ == "__main__":
    main()
