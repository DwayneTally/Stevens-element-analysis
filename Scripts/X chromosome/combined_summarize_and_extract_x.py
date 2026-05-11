#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path
import pandas as pd


def read_list_file(path):
    if path is None:
        return set()
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def find_file(genome, directory, exts):
    for ext in exts:
        p = directory / f"{genome}{ext}"
        if p.exists():
            return p
    return None


def read_fasta_ids_and_subset(fasta_file, keep_ids, out_file):
    n_written = 0
    with open_text(fasta_file) as fin, open(out_file, "w") as fout:
        write_block = False
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                write_block = header in keep_ids
                if write_block:
                    fout.write(line)
                    n_written += 1
            else:
                if write_block:
                    fout.write(line)
    return n_written


def main():
    parser = argparse.ArgumentParser(
        description="Run summarize_x_genes and then subset GENESPACE bed/peptide files to inferred X chromosomes"
    )

    #summarize_x_genes.py arguments
    parser.add_argument("--bed", required=True)
    parser.add_argument("--pangenes", required=True)
    parser.add_argument("--x-chr", default="NC_007416.3")
    parser.add_argument("--split-genomes", default=None)
    parser.add_argument("--exclude-genomes", default=None)
    parser.add_argument("--out-prefix", default="tribolium_X")

    #genespace_extract_x.py arguments
    parser.add_argument("--genespace-root", required=True, help="GENESPACE run directory containing bed/ and peptide/")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--genome-col", default="genome", help="Genome column in summary CSV")
    parser.add_argument("--x-col", default="x_chromosome", help="X chromosome column in summary CSV")

    args = parser.parse_args()

    split_genomes = read_list_file(args.split_genomes)
    exclude_genomes = read_list_file(args.exclude_genomes)

    bed = pd.read_csv(args.bed, sep="\t", header=None, dtype=str)

    trib_x_genes = set(bed.loc[bed[0] == args.x_chr, 3])

    print("Tribolium X genes:", len(trib_x_genes))

    pg = pd.read_csv(args.pangenes, sep="\t", compression="gzip", dtype=str)

    trib_x_pangenes = set(pg.loc[pg["id"].isin(trib_x_genes), "pgID"])

    print("Tribolium X pangenes:", len(trib_x_pangenes))

    hits = pg.loc[pg["pgID"].isin(trib_x_pangenes), ["genome", "chr", "id", "pgID"]].copy()
    hits.columns = ["genome", "chr", "gene_id", "pgID"]

    if exclude_genomes:
        hits = hits.loc[~hits["genome"].isin(exclude_genomes)]

    hits.to_csv(f"{args.out_prefix}_pangene_hits.csv", index=False)

    chr_counts = (
        hits.groupby(["genome", "chr"])
        .size()
        .reset_index(name="chromosome_x_genes")
        .sort_values(["genome", "chromosome_x_genes"], ascending=[True, False])
    )

    chr_counts.to_csv(f"{args.out_prefix}_chr_counts.csv", index=False)

    totals = (
        hits.groupby("genome")["pgID"]
        .nunique()
        .reset_index(name="number_pangenes_found")
    )

    summary_rows = []

    split_present = False

    for genome, sub in chr_counts.groupby("genome"):

        sub = sub.sort_values("chromosome_x_genes", ascending=False).reset_index(drop=True)

        top_chr = sub.loc[0, "chr"]
        top_count = int(sub.loc[0, "chromosome_x_genes"])

        second_chr = ""
        second_count = ""

        if genome in split_genomes and len(sub) > 1:
            split_present = True
            second_chr = sub.loc[1, "chr"]
            second_count = int(sub.loc[1, "chromosome_x_genes"])
            x_chr = f"{top_chr};{second_chr}"
        else:
            x_chr = top_chr

        total_pg = int(
            totals.loc[totals["genome"] == genome, "number_pangenes_found"].iloc[0]
        )

        summary_rows.append(
            {
                "genome": genome,
                "x_chromosome": x_chr,
                "chromosome_x_genes1": top_count,
                "chromosome_x_genes2": second_count,
                "number_pangenes_found": total_pg,
            }
        )

    summary = pd.DataFrame(summary_rows)

    if split_present:

        summary = summary[
            [
                "genome",
                "x_chromosome",
                "chromosome_x_genes1",
                "chromosome_x_genes2",
                "number_pangenes_found",
            ]
        ]

    else:

        summary["chromosome_x_genes"] = summary["chromosome_x_genes1"]

        summary = summary[
            [
                "genome",
                "x_chromosome",
                "chromosome_x_genes",
                "number_pangenes_found",
            ]
        ]

    summary_file = f"{args.out_prefix}_homology_summary.csv"
    summary.to_csv(summary_file, index=False)

    print("\nOutput:")
    print(summary.to_string(index=False))

    root = Path(args.genespace_root)
    bed_dir = root / "bed"
    pep_dir = root / "peptide"
    outdir = Path(args.outdir)
    out_bed = outdir / "bed"
    out_pep = outdir / "peptide"
    outdir.mkdir(parents=True, exist_ok=True)
    out_bed.mkdir(exist_ok=True)
    out_pep.mkdir(exist_ok=True)

    bed_exts = [".bed"]
    pep_exts = [".fa", ".faa", ".fasta", ".fa.gz", ".faa.gz", ".fasta.gz"]

    subset_summary_rows = []

    with open(summary_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            genome = row[args.genome_col].strip()
            x_chrom = row[args.x_col].strip()

            if not genome or not x_chrom:
                continue

            x_chroms = {x.strip() for x in x_chrom.split(";") if x.strip()}

            bed_file = find_file(genome, bed_dir, bed_exts)
            pep_file = find_file(genome, pep_dir, pep_exts)

            status = []
            kept_gene_ids = set()
            bed_out_file = out_bed / f"{genome}.bed"
            pep_out_file = out_pep / f"{genome}.fa"

            if bed_file is None:
                status.append("missing_bed")
            else:
                n_bed_rows = 0
                with open_text(bed_file) as fin, open(bed_out_file, "w") as fout:
                    for line in fin:
                        if not line.strip():
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 4:
                            continue
                        seqid = parts[0]
                        gene_id = parts[3]
                        if seqid in x_chroms:
                            fout.write(line)
                            kept_gene_ids.add(gene_id)
                            n_bed_rows += 1
                status.append("ok_bed")

            if pep_file is None:
                status.append("missing_peptide")
                n_fasta = 0
            else:
                n_fasta = read_fasta_ids_and_subset(pep_file, kept_gene_ids, pep_out_file)
                status.append("ok_peptide")

            subset_summary_rows.append({
                "genome": genome,
                "x_chromosome": x_chrom,
                "bed_file": str(bed_file) if bed_file else "",
                "peptide_file": str(pep_file) if pep_file else "",
                "n_x_bed_rows": len(kept_gene_ids),
                "n_x_peptides": n_fasta,
                "status": ";".join(status),
            })

    subset_summary_out = outdir / "subset_summary.csv"
    with open(subset_summary_out, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "genome", "x_chromosome", "bed_file", "peptide_file",
                "n_x_bed_rows", "n_x_peptides", "status"
            ],
        )
        writer.writeheader()
        writer.writerows(subset_summary_rows)

    print(f"\nWrote X-only BED files to: {out_bed}")
    print(f"Wrote X-only peptide FASTAs to: {out_pep}")
    print(f"Wrote summary to: {subset_summary_out}")


if __name__ == "__main__":
    main()
