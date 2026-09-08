#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import defaultdict

ELEMENTS = ["A", "C", "E", "G", "H", "X"]

def read_fasta(path):
    """
    Simple FASTA parser.
    Yields:
        header_id, full_header, sequence
    """

    header = None
    seq_parts = []

    with open(path) as f:

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if header is not None:

                    full_header = header
                    header_id = full_header.split()[0]

                    yield (
                        header_id,
                        full_header,
                        "".join(seq_parts)
                    )

                header = line[1:]
                seq_parts = []

            else:

                seq_parts.append(
                    line.strip()
                )

        if header is not None:

            full_header = header
            header_id = full_header.split()[0]

            yield (
                header_id,
                full_header,
                "".join(seq_parts)
            )



def write_fasta_record(out, header, seq, width=80):

    out.write(f">{header}\n")

    for i in range(
        0,
        len(seq),
        width
    ):

        out.write(
            seq[i:i + width] + "\n"
        )


# READ CHROMOSOME, STEVENS ASSIGNMENTS

def read_assignments(path, elements):

    """
    Expected columns from gene_movement_busco_v2.py:

    Species
    Chromosome
    Chromosome_element
    ...
    Assignment_flag

    Returns:
        element -> species -> set(chromosomes)
    """

    mapping = defaultdict(
        lambda: defaultdict(set)
    )

    with open(path) as f:

        header = f.readline().rstrip("\n").split("\t")

        col = {
            name: i
            for i, name in enumerate(header)
        }

        required = [
            "Species",
            "Chromosome",
            "Chromosome_element"
        ]

        for r in required:

            if r not in col:

                raise ValueError(
                    f"Missing required column in "
                    f"assignment file: {r}"
                )

        for line in f:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            species = fields[
                col["Species"]
            ].strip()

            chromosome = fields[
                col["Chromosome"]
            ].strip()

            element = fields[
                col["Chromosome_element"]
            ].strip()

            if element not in elements:
                continue

            mapping[element][species].add(
                chromosome
            )

    return mapping


# READ BED AND GET GENES ON CHROMOSOMES

def subset_bed(
    bed_path,
    output_path,
    target_chromosomes
):

    gene_ids = set()

    n_total = 0
    n_kept = 0

    with open(bed_path) as infile, \
         open(output_path, "w") as outfile:

        for line in infile:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 4:
                continue

            n_total += 1

            chromosome = fields[0]
            gene_id = fields[3]

            if chromosome not in target_chromosomes:
                continue

            outfile.write(
                line + "\n"
            )

            gene_ids.add(
                gene_id
            )

            n_kept += 1

    return gene_ids, n_total, n_kept

def subset_peptides(
    fasta_path,
    output_path,
    gene_ids
):

    fasta_ids = set()

    n_total = 0
    n_kept = 0

    with open(output_path, "w") as outfile:

        for header_id, full_header, seq in read_fasta(
            fasta_path
        ):

            n_total += 1

            if header_id not in gene_ids:
                continue

            write_fasta_record(
                outfile,
                full_header,
                seq
            )

            fasta_ids.add(
                header_id
            )

            n_kept += 1

    missing_ids = (
        gene_ids - fasta_ids
    )

    return (
        n_total,
        n_kept,
        missing_ids
    )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate Stevens-element-specific BED and peptide "
            "datasets for GENESPACE using chromosome assignments."
        )
    )

    parser.add_argument(
        "--assignments",
        required=True,
        help=(
            "chromosome_stevens_assignments.tsv from "
            "gene_movement_busco_v2.py"
        )
    )

    parser.add_argument(
        "--bed-dir",
        required=True,
        help="Master BED directory"
    )

    parser.add_argument(
        "--peptide-dir",
        required=True,
        help="Master peptide FASTA directory"
    )

    parser.add_argument(
        "--outdir",
        default="genespace_by_element"
    )

    parser.add_argument(
        "--elements",
        nargs="+",
        default=ELEMENTS,
        help=(
            "Stevens elements to generate. "
            "Default: A C E G H X"
        )
    )

    args = parser.parse_args()

    bed_dir = Path(
        args.bed_dir
    )

    peptide_dir = Path(
        args.peptide_dir
    )

    outdir = Path(
        args.outdir
    )

    elements = set(
        args.elements
    )

    # READ CHROMOSOME ASSIGNMENTS

    assignment_map = read_assignments(
        args.assignments,
        elements
    )

    print()
    print(
        "GENESPACE ELEMENT-SPECIFIC DATASET"
    )
 
    print()
    print(
        "Elements:"
    )

    print(
        "  " + " ".join(
            sorted(elements)
        )
    )

    bed_files = {
        p.stem: p
        for p in bed_dir.glob("*.bed")
    }

    peptide_files = {
        p.stem: p
        for p in peptide_dir.glob("*.fa")
    }

    shared_species = sorted(
        set(bed_files)
        &
        set(peptide_files)
    )

    print()
    print(
        f"Species with matching BED + peptide: "
        f"{len(shared_species)}"
    )

    if not shared_species:

        raise RuntimeError(
            "No matching BED/peptide filenames found."
        )

    report_path = (
        outdir /
        "element_dataset_report.tsv"
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(report_path, "w") as report:

        report.write(
            "Element\t"
            "Species\t"
            "Chromosomes\t"
            "BED_genes_total\t"
            "BED_genes_kept\t"
            "Peptides_total\t"
            "Peptides_kept\t"
            "Missing_peptides\n"
        )

        # PROCESS EACH ELEMENT

        for element in sorted(
            elements
        ):

            print()
            print(
                f"STEVENS {element}"
            )

            element_dir = (
                outdir /
                element
            )

            element_bed_dir = (
                element_dir /
                "bed"
            )

            element_peptide_dir = (
                element_dir /
                "peptide"
            )

            element_bed_dir.mkdir(
                parents=True,
                exist_ok=True
            )

            element_peptide_dir.mkdir(
                parents=True,
                exist_ok=True
            )

            n_species_written = 0

            # PROCESS EACH SPECIES


            for base in shared_species:
                # BED/FASTA basename contains accession.

                matching_species = None

                for species_name in (
                    assignment_map[
                        element
                    ].keys()
                ):

                    if base.startswith(
                        species_name + "_GCA_"
                    ) or base.startswith(
                        species_name + "_GCF_"
                    ):

                        matching_species = (
                            species_name
                        )

                        break

                    if base == species_name:

                        matching_species = (
                            species_name
                        )

                        break

                if matching_species is None:
                    continue

                chromosomes = (
                    assignment_map[
                        element
                    ][
                        matching_species
                    ]
                )

                if not chromosomes:
                    continue

                bed_input = (
                    bed_files[
                        base
                    ]
                )

                pep_input = (
                    peptide_files[
                        base
                    ]
                )

                bed_output = (
                    element_bed_dir /
                    f"{base}.bed"
                )

                pep_output = (
                    element_peptide_dir /
                    f"{base}.fa"
                )

                (
                    gene_ids,
                    bed_total,
                    bed_kept
                ) = subset_bed(
                    bed_input,
                    bed_output,
                    chromosomes
                )

                (
                    pep_total,
                    pep_kept,
                    missing_peptides
                ) = subset_peptides(
                    pep_input,
                    pep_output,
                    gene_ids
                )

                if bed_kept == 0:

                    print(
                        f"WARNING: {base} {element}: "
                        f"0 BED genes found on "
                        f"{','.join(sorted(chromosomes))}"
                    )

                if missing_peptides:

                    print(
                        f"WARNING: {base} {element}: "
                        f"{len(missing_peptides)} BED gene IDs "
                        f"missing from peptide FASTA"
                    )

                report.write(
                    f"{element}\t"
                    f"{base}\t"
                    f"{','.join(sorted(chromosomes))}\t"
                    f"{bed_total}\t"
                    f"{bed_kept}\t"
                    f"{pep_total}\t"
                    f"{pep_kept}\t"
                    f"{len(missing_peptides)}\n"
                )

                print(
                    f"{base}: "
                    f"{bed_kept} genes "
                    f"from {len(chromosomes)} chromosome(s)"
                )

                n_species_written += 1

            print()
            print(
                f"Species written for {element}: "
                f"{n_species_written}"
            )


    print()
    print(
        f"Output directory:"
    )

    print(
        f"  {outdir}"
    )

    print()
    print(
        f"Report:"
    )

    print(
        f"  {report_path}"
    )


if __name__ == "__main__":
    main()
