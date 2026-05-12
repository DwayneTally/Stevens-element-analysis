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

        gffs = list(summary_dir.glob("*.filteredRepeats.gff"))
        fastas = list(summary_dir.glob("*.softmasked.fasta"))

        if gffs and fastas:
            mapping[norm_name(genome_dir.name)] = {
                "gff": gffs[0],
                "fasta": fastas[0],
            }

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


def get_fasta_lengths(fasta_file, genome):
    lengths = {}
    current = None

    with open(fasta_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                lengths[current] = 0
                continue

            if current is not None:
                seq = line.strip()
                lengths[current] += sum(1 for ch in seq if ch.isalpha())

    # species-specific cleanup
    if genome == "Crioceris_asparagi_GCA_958507055.1":
        lengths = {
            k: v for k, v in lengths.items()
            if k.startswith("OY293") and k != "OY293834.1"
        }

    return lengths


def collapse_label(label):
    if label.startswith("DNA"):
        return "DNA"
    if label.startswith("LINE"):
        return "LINE"
    if label.startswith("LTR"):
        return "LTR"
    if label.startswith("RC"):
        return "RC"
    if label.startswith("SINE"):
        return "SINE"
    if label == "Satellite":
        return "Satellite"
    if label == "Simple_repeat":
        return "Simple_repeat"
    if label == "Low_complexity":
        return "Low_complexity"
    if label == "Unknown":
        return "Unknown"
    return "Unknown"


def parse_gff_bp_by_group(gff_file, x_scaffolds, allowed_scaffolds):
    categories = [
        "DNA", "LINE", "LTR", "RC", "SINE",
        "Satellite", "Simple_repeat", "Low_complexity", "Unknown"
    ]

    counts = {
        "X_noNeoX": {c: 0 for c in categories},
        "Autosomes": {c: 0 for c in categories},
    }

    with open(gff_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            seqid = parts[0]
            label = parts[2]
            start = int(parts[3])
            end = int(parts[4])

            if seqid not in allowed_scaffolds:
                continue

            group = "X_noNeoX" if seqid in x_scaffolds else "Autosomes"
            cat = collapse_label(label)
            bp = end - start + 1
            counts[group][cat] += bp

    return counts


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse EarlGrey filteredRepeats.gff and compute X vs autosome repeat "
            "statistics, including TE proportion of total repeats."
        )
    )
    parser.add_argument("--earlgrey-base", required=True,
                        help="Directory containing one EarlGrey genome directory per species")
    parser.add_argument("--x-bed-dir", required=True,
                        help="Directory containing BED files of ancestral X scaffolds")
    parser.add_argument("--exclude-list", default=None,
                        help="Optional text file listing genomes to exclude")
    parser.add_argument("--out", default="x_vs_autosomes_repeat_categories.csv",
                        help="Output CSV filename")
    args = parser.parse_args()

    earlgrey_map = build_earlgrey_map(Path(args.earlgrey_base))
    bed_map = build_bed_map(Path(args.x_bed_dir))
    exclude = read_exclude_list(args.exclude_list)

    categories = [
        "DNA", "LINE", "LTR", "RC", "SINE",
        "Satellite", "Simple_repeat", "Low_complexity", "Unknown"
    ]

    rows = []
    missing_earlgrey = []
    skipped_missing_x = []
    excluded = []

    for genome in sorted(bed_map.keys()):
        if genome in exclude:
            excluded.append(genome)
            continue

        if genome not in earlgrey_map:
            missing_earlgrey.append(genome)
            continue

        gff_file = earlgrey_map[genome]["gff"]
        fasta_file = earlgrey_map[genome]["fasta"]
        bed_file = bed_map[genome]

        x_scaffolds = read_x_scaffolds_from_bed(bed_file)
        fasta_lengths = get_fasta_lengths(fasta_file, genome)
        allowed_scaffolds = set(fasta_lengths.keys())

        matched_x = x_scaffolds & allowed_scaffolds
        if len(matched_x) == 0:
            skipped_missing_x.append(genome)
            print(f"WARNING: skipping {genome}; no X scaffolds found in FASTA")
            continue

        x_len = sum(fasta_lengths[s] for s in allowed_scaffolds if s in x_scaffolds)
        auto_len = sum(fasta_lengths[s] for s in allowed_scaffolds if s not in x_scaffolds)

        if x_len == 0 or auto_len == 0:
            print(f"WARNING: skipping {genome}; zero X or autosome length")
            continue

        bp_counts = parse_gff_bp_by_group(gff_file, x_scaffolds, allowed_scaffolds)

        for group, comp_len in [("X_noNeoX", x_len), ("Autosomes", auto_len)]:
            row = {
                "genome": genome,
                "group": group,
                "compartment_bp": comp_len,
            }

            total_masked = 0
            for cat in categories:
                bp = bp_counts[group][cat]
                pct = round((bp / comp_len) * 100, 2) if comp_len > 0 else 0.0
                row[f"{cat}_bp"] = bp
                row[f"{cat}_percent"] = pct
                total_masked += bp

            row["masked_bp"] = total_masked
            row["masked_percent"] = round((total_masked / comp_len) * 100, 2) if comp_len > 0 else 0.0

            te_bp = (
                row["DNA_bp"] +
                row["LINE_bp"] +
                row["LTR_bp"] +
                row["RC_bp"] +
                row["SINE_bp"] +
                row["Unknown_bp"]
            )
            row["te_bp"] = te_bp

            # Percent of compartment occupied by TE-like repeats
            row["te_percent_of_compartment"] = round((te_bp / comp_len) * 100, 2) if comp_len > 0 else 0.0

            # Percent of total repeat content that is TE-like
            row["te_percent_of_repeats"] = round((te_bp / total_masked) * 100, 2) if total_masked > 0 else 0.0

            rows.append(row)

    fieldnames = [
        "genome", "group", "compartment_bp",
        "DNA_bp", "DNA_percent",
        "LINE_bp", "LINE_percent",
        "LTR_bp", "LTR_percent",
        "RC_bp", "RC_percent",
        "SINE_bp", "SINE_percent",
        "Satellite_bp", "Satellite_percent",
        "Simple_repeat_bp", "Simple_repeat_percent",
        "Low_complexity_bp", "Low_complexity_percent",
        "Unknown_bp", "Unknown_percent",
        "masked_bp", "masked_percent",
        "te_bp", "te_percent_of_compartment", "te_percent_of_repeats"
    ]

    with open(args.out, "w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.out}")

    if missing_earlgrey:
        print("\nMissing EarlGrey files:")
        for g in missing_earlgrey:
            print(f"  {g}")

    if skipped_missing_x:
        print("\nSkipped because X scaffolds missing from FASTA:")
        for g in skipped_missing_x:
            print(f"  {g}")

    if excluded:
        print("\nExcluded by exclude-list:")
        for g in excluded:
            print(f"  {g}")


if __name__ == "__main__":
    main()
