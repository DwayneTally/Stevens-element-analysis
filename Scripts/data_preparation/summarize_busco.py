#!/usr/bin/env python3
import os
import glob
import re
import csv
import sys

def parse_short_summary(path):
    #Parse BUSCO short_summary*.txt file and return a dict with key metrics.
    lineage = None
    mode = None
    busco_version = None
    stats_line = None

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("# The lineage dataset is:"):
                # e.g. "# The lineage dataset is: endopterygota_odb10 (Creation date: ..., number of BUSCOs: 255)"
                lineage = line.split("is:")[1].split("(")[0].strip()
            elif line.startswith("# BUSCO was run in mode:"):
                mode = line.split("mode:")[1].strip()
            elif line.startswith("# BUSCO version is:"):
                busco_version = line.split("is:")[1].strip()
            elif line.startswith("C:"):
                stats_line = line

    if stats_line is None:
        raise ValueError(f"No 'C:' stats line found in {path}")

    #Example:
    #C:98.7%[S:97.5%,D:1.2%],F:0.5%,M:0.8%,n:1367
    pattern = re.compile(
        r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],"
        r"F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)"
    )
    m = pattern.search(stats_line)
    if not m:
        raise ValueError(f"Could not parse stats line in {path}: {stats_line}")

    return {
        "lineage": lineage,
        "mode": mode,
        "busco_version": busco_version,
        "C": float(m.group("C")),
        "S": float(m.group("S")),
        "D": float(m.group("D")),
        "F": float(m.group("F")),
        "M": float(m.group("M")),
        "n": int(m.group("n")),
    }

def parse_ids_from_busco_dir(busco_dir_name):
    """
    return:
      genome_id = Abax_parallelepipedus_GCA_964197645.1_genomic
      species   = "Abax parallelepipedus"
      accession = "GCA_964197645.1"
    """
    genome_id = busco_dir_name
    if genome_id.startswith("busco_"):
        genome_id = genome_id[len("busco_"):]

    m = re.search(r"(GC[AF]_\d+\.\d+)", genome_id)
    accession = m.group(1) if m else ""

    species = ""
    if accession:
        species_part = genome_id.split("_" + accession)[0]
        species = species_part.replace("_", " ")
    else:
        # fallback: take first two tokens as species-like name
        toks = genome_id.split("_")
        if len(toks) >= 2:
            species = " ".join(toks[:2])
        else:
            species = genome_id

    return genome_id, species, accession

def main(root_dir, outfile):
    pattern = os.path.join(root_dir, "busco_*", "short_summary.specific.*.txt")
    paths = sorted(glob.glob(pattern))

    if not paths:
        sys.stderr.write(f"No short_summary.specific.*.txt files found under {root_dir}\n")
        sys.exit(1)

    fieldnames = [
        "genome_id",
        "species",
        "assembly_accession",
        "busco_dir",
        "lineage",
        "mode",
        "busco_version",
        "n",
        "C",
        "S",
        "D",
        "F",
        "M",
    ]

    with open(outfile, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for path in paths:
            busco_dir = os.path.basename(os.path.dirname(path))
            genome_id, species, accession = parse_ids_from_busco_dir(busco_dir)
            stats = parse_short_summary(path)

            row = {
                "genome_id": genome_id,
                "species": species,
                "assembly_accession": accession,
                "busco_dir": busco_dir,
                "lineage": stats["lineage"],
                "mode": stats["mode"],
                "busco_version": stats["busco_version"],
                "n": stats["n"],
                "C": stats["C"],
                "S": stats["S"],
                "D": stats["D"],
                "F": stats["F"],
                "M": stats["M"],
            }
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(
            f"Usage: {sys.argv[0]} /path/to/busco_beetle_result busco_beetle_summary.tsv\n"
        )
        sys.exit(1)

    root = sys.argv[1]
    out = sys.argv[2]
    main(root, out)

