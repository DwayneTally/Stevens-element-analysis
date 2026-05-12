#!/usr/bin/env python3

import os
from pathlib import Path

def extract_pep_ids(peptide_fasta):
    ids = set()
    with open(peptide_fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids

def parse_gff3_to_bed(gff3_path, pep_ids, output_bed):
    with open(gff3_path) as gff, open(output_bed, "w") as bed:
        for line in gff:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) != 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attributes = cols
            if feature != "mRNA":
                continue
            attr_dict = dict(field.split("=", 1) for field in attributes.split(";") if "=" in field)
            transcript_id = attr_dict.get("ID")
            if transcript_id and transcript_id in pep_ids:
                bed.write(f"{chrom}\t{int(start)-1}\t{end}\t{transcript_id}\n")

def main():
    gff_dir = Path("genomes")
    pep_dir = Path("peptide")
    output_dir = Path("bed")
    output_dir.mkdir(parents=True, exist_ok=True)

    gff_files = sorted(gff_dir.glob("*.gff"))

    for gff_file in gff_files:
        base = gff_file.stem  # e.g., "Tribolium_castaneum_GCF_000002335.3"
        pep_file = pep_dir / f"{base}.fa"
        bed_file = output_dir / f"{base}.bed"

        if not pep_file.exists():
            print(f"Skipping {base} — peptide file not found.")
            continue

        pep_ids = extract_pep_ids(pep_file)
        parse_gff3_to_bed(gff_file, pep_ids, bed_file)
        print(f"Wrote: {bed_file}")

if __name__ == "__main__":
    main()

