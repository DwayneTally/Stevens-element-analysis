#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

TARGET_CLASSES = {
    "DNA": "DNA",
    "Rolling Circle": "Rolling Circle",
    "Penelope": "Penelope",
    "LINE": "LINE",
    "SINE": "SINE",
    "LTR": "LTR",
    "Other (Simple Repeat, Microsatellite, RNA)": "Simple Repeat (Other)",
    "Unclassified": "Unclassified",
}

NON_REPEAT_LABEL = "Non-Repeat"


def parse_highlevelcount(file_path):
    """
    Parse one EarlGrey *.highLevelCount.txt file.

    Returns:
      genome_size (int or None)
      non_repeat_coverage (int or None)
      class_coverages (dict)
    """
    genome_size = None
    non_repeat_coverage = None
    class_coverages = {v: 0 for v in TARGET_CLASSES.values()}

    with open(file_path, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    if not lines:
        return genome_size, non_repeat_coverage, class_coverages

    for line in lines[1:]:  # skip header
        parts = re.split(r"\t+|\s{2,}", line.strip())

        if len(parts) < 6:
            continue

        te_class = parts[0].strip()
        coverage_bp = parts[1].strip()
        genome_size_field = parts[4].strip()

        try:
            coverage_bp = int(float(coverage_bp))
        except ValueError:
            continue

        try:
            genome_size_val = int(float(genome_size_field))
        except ValueError:
            genome_size_val = None

        if genome_size is None and genome_size_val is not None:
            genome_size = genome_size_val

        if te_class == NON_REPEAT_LABEL:
            non_repeat_coverage = coverage_bp
            continue

        # skip nested rows
        if te_class.endswith("-nested"):
            continue

        if te_class in TARGET_CLASSES:
            out_name = TARGET_CLASSES[te_class]
            class_coverages[out_name] = coverage_bp

    return genome_size, non_repeat_coverage, class_coverages


def find_highlevelcount_files(base_dir):
    files = []
    for genome_dir in sorted(base_dir.iterdir()):
        if not genome_dir.is_dir():
            continue

        summary_dir = genome_dir / "summaryFiles"
        if not summary_dir.exists():
            continue

        matches = list(summary_dir.glob("*.highLevelCount.txt"))
        if matches:
            files.extend(matches)

    return files


def genome_name_from_path(file_path):
    return file_path.parent.parent.name


def main():
    parser = argparse.ArgumentParser(
        description="Parse EarlGrey highLevelCount.txt files and compute TE composition percentages plus whole-genome masked/unmasked percentages."
    )
    parser.add_argument(
        "--base-dir",
        required=True,
        help="Directory containing genome subdirectories with summaryFiles/"
    )
    parser.add_argument(
        "--out",
        default="earlgrey_repeat_percentages.csv",
        help="Output CSV file"
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    out_csv = Path(args.out)

    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")

    files = find_highlevelcount_files(base_dir)
    if not files:
        raise FileNotFoundError(f"No *.highLevelCount.txt files found under {base_dir}")

    output_rows = []

    for file_path in files:
        genome = genome_name_from_path(file_path)
        genome_size, non_repeat_coverage, class_coverages = parse_highlevelcount(file_path)

        if genome_size is None:
            print(f"WARNING: could not determine Genome Size for {genome}, skipping")
            continue

        if non_repeat_coverage is None:
            print(f"WARNING: could not determine Non-Repeat coverage for {genome}, skipping")
            continue

        total_te_bp = sum(class_coverages.values())

        if total_te_bp <= 0:
            print(f"WARNING: total selected TE bp <= 0 for {genome}, skipping")
            continue

        masked_percent = round(((genome_size - non_repeat_coverage) / genome_size) * 100, 2)
        unmasked_percent = round((non_repeat_coverage / genome_size) * 100, 2)

        row = {
            "Genome": genome,
            "Genome Size": genome_size,
            "DNA": round((class_coverages["DNA"] / total_te_bp) * 100, 2),
            "Rolling Circle": round((class_coverages["Rolling Circle"] / total_te_bp) * 100, 2),
            "Penelope": round((class_coverages["Penelope"] / total_te_bp) * 100, 2),
            "LINE": round((class_coverages["LINE"] / total_te_bp) * 100, 2),
            "SINE": round((class_coverages["SINE"] / total_te_bp) * 100, 2),
            "LTR": round((class_coverages["LTR"] / total_te_bp) * 100, 2),
            "Simple Repeat (Other)": round((class_coverages["Simple Repeat (Other)"] / total_te_bp) * 100, 2),
            "Unclassified": round((class_coverages["Unclassified"] / total_te_bp) * 100, 2),
            "Masked": masked_percent,
            "Unmasked": unmasked_percent,
        }

        output_rows.append(row)

    fieldnames = [
        "Genome",
        "Genome Size",
        "DNA",
        "Rolling Circle",
        "Penelope",
        "LINE",
        "SINE",
        "LTR",
        "Simple Repeat (Other)",
        "Unclassified",
        "Masked",
        "Unmasked",
    ]

    with open(out_csv, "w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} genomes to {out_csv}")


if __name__ == "__main__":
    main()
