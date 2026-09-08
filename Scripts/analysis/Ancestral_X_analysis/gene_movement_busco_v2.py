#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


FOCAL_ELEMENTS = ["X", "E", "A", "G", "C", "H"]

def read_stevens(path):

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str
    )

    df.columns = [
        x.strip()
        for x in df.columns
    ]

    if "Orthogroups" not in df.columns:
        raise ValueError(
            "stevens.tsv must contain an 'Orthogroups' column"
        )

    if "Stevens" not in df.columns:
        raise ValueError(
            "stevens.tsv must contain a 'Stevens' column"
        )

    df = df[
        ["Orthogroups", "Stevens"]
    ].copy()

    df.columns = [
        "BUSCO",
        "Ancestral_element"
    ]

    df["BUSCO"] = (
        df["BUSCO"]
        .astype(str)
        .str.strip()
    )

    df["Ancestral_element"] = (
        df["Ancestral_element"]
        .astype(str)
        .str.strip()
    )

    df = df.drop_duplicates()

    # Identify BUSCOs assigned to >1 Stevens element

    assignment_counts = (
        df.groupby("BUSCO")[
            "Ancestral_element"
        ]
        .nunique()
    )

    conflicting_buscos = set(
        assignment_counts[
            assignment_counts > 1
        ].index
    )

    if conflicting_buscos:

        print()
        print(
            "WARNING: CONFLICTING STEVENS ASSIGNMENTS"
        )

        for busco in sorted(
            conflicting_buscos
        ):

            elements = sorted(
                df.loc[
                    df["BUSCO"] == busco,
                    "Ancestral_element"
                ].unique()
            )

            print(
                f"  {busco}: "
                f"{','.join(elements)}"
            )

        print()
        print(
            f"Excluding "
            f"{len(conflicting_buscos)} "
            f"ambiguous BUSCO IDs."
        )
        print()

        df = df[
            ~df["BUSCO"].isin(
                conflicting_buscos
            )
        ].copy()

    df = df.drop_duplicates(
        subset=["BUSCO"]
    )

    print(
        f"Usable Stevens-mapped BUSCOs: "
        f"{len(df)}"
    )

    return df

def read_busco_table(path):

    rows = []

    with open(path) as f:

        for line in f:

            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                continue

            busco = fields[0].strip()
            status = fields[1].strip()
            sequence = fields[2].strip()

            rows.append(
                {
                    "BUSCO": busco,
                    "Status": status,
                    "Chromosome": sequence,
                }
            )

    df = pd.DataFrame(rows)

    if df.empty:
        return df

    # Keep only single-copy Complete BUSCOs
    df = df[
        df["Status"] == "Complete"
    ].copy()

    # Exclude accidental repeated Complete BUSCO rows

    duplicated = df[
        df["BUSCO"].duplicated(
            keep=False
        )
    ]

    if len(duplicated) > 0:

        bad_ids = set(
            duplicated["BUSCO"]
        )

        print(
            f"WARNING: found "
            f"{len(bad_ids)} BUSCO IDs with "
            f">1 Complete row; excluding them."
        )

        df = df[
            ~df["BUSCO"].isin(
                bad_ids
            )
        ].copy()

    return df


# ASSIGN CHROMOSOMES TO STEVENS ELEMENTS

def assign_chromosomes(dat, species):

    """
    Assign each chromosome to the Stevens element represented
    by the largest number of BUSCOs.

    IMPORTANT:
    Multiple chromosomes are allowed to receive the same
    Stevens identity. This accommodates chromosome fissions.

    Also records the second-most-common Stevens element and
    the difference between first and second place.
    """

    counts = (
        dat
        .groupby(
            [
                "Chromosome",
                "Ancestral_element"
            ]
        )
        .size()
        .reset_index(
            name="N"
        )
    )

    output_rows = []

    for chromosome, sub in counts.groupby(
        "Chromosome"
    ):

        sub = (
            sub
            .sort_values(
                "N",
                ascending=False
            )
            .reset_index(
                drop=True
            )
        )

        total = int(
            sub["N"].sum()
        )

        top_element = (
            sub.loc[
                0,
                "Ancestral_element"
            ]
        )

        top_n = int(
            sub.loc[
                0,
                "N"
            ]
        )

        if len(sub) >= 2:

            second_element = (
                sub.loc[
                    1,
                    "Ancestral_element"
                ]
            )

            second_n = int(
                sub.loc[
                    1,
                    "N"
                ]
            )

        else:

            second_element = "NA"
            second_n = 0

        top_proportion = (
            top_n / total
            if total > 0
            else 0
        )

        second_proportion = (
            second_n / total
            if total > 0
            else 0
        )

        dominance_margin = (
            top_proportion
            -
            second_proportion
        )

        if second_n > 0:

            top_to_second_ratio = (
                top_n / second_n
            )

        else:

            top_to_second_ratio = float("inf")

        output_rows.append(
            {
                "Species":
                    species,

                "Chromosome":
                    chromosome,

                "Chromosome_element":
                    top_element,

                "Dominant_BUSCOs":
                    top_n,

                "Total_BUSCOs":
                    total,

                "Dominance":
                    top_proportion,

                "Second_element":
                    second_element,

                "Second_BUSCOs":
                    second_n,

                "Second_proportion":
                    second_proportion,

                "Dominance_margin":
                    dominance_margin,

                "Top_to_second_ratio":
                    top_to_second_ratio,
            }
        )

    return pd.DataFrame(
        output_rows
    )

def main():

    parser = argparse.ArgumentParser(
        description=(
            "BUSCO-based analysis of gene retention and "
            "movement among Stevens elements while allowing "
            "multiple descendant chromosomes per ancestral "
            "element."
        )
    )

    parser.add_argument(
        "--species",
        required=True,
        help=(
            "TSV with columns: Species and full_table"
        )
    )

    parser.add_argument(
        "--stevens",
        required=True,
        help=(
            "Stevens BUSCO mapping used by Element_vis"
        )
    )

    parser.add_argument(
        "--outdir",
        default="gene_movement_results_v2"
    )

    parser.add_argument(
        "--tie-margin",
        type=float,
        default=0.10,
        help=(
            "Flag chromosomes where the proportional "
            "difference between the top and second Stevens "
            "element is < this value. "
            "Default = 0.10. "
            "These chromosomes are still analyzed."
        )
    )

    args = parser.parse_args()

    outdir = Path(
        args.outdir
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    stevens = read_stevens(
        args.stevens
    )

    species_table = pd.read_csv(
        args.species,
        sep="\t",
        dtype=str
    )

    required_columns = {
        "Species",
        "full_table"
    }

    if not required_columns.issubset(
        species_table.columns
    ):

        raise ValueError(
            "Species TSV must contain columns:\n"
            "Species\tfull_table"
        )

    all_gene_rows = []
    chromosome_rows = []
    composition_rows = []

    for _, row in species_table.iterrows():

        species = row[
            "Species"
        ].strip()

        full_table = Path(
            row[
                "full_table"
            ].strip()
        )

        print()
        print(
            species
        )

        if not full_table.exists():

            print(
                f"WARNING: missing BUSCO table:"
            )

            print(
                f"  {full_table}"
            )

            continue


        busco = read_busco_table(
            full_table
        )

        print(
            f"Complete single-copy BUSCOs: "
            f"{len(busco)}"
        )

        dat = busco.merge(
            stevens,
            on="BUSCO",
            how="inner"
        )

        print(
            f"BUSCOs with usable Stevens assignment: "
            f"{len(dat)}"
        )

        if dat.empty:

            print(
                "WARNING: no BUSCOs matched Stevens mapping."
            )

            continue

        composition = (
            dat
            .groupby(
                [
                    "Chromosome",
                    "Ancestral_element"
                ]
            )
            .size()
            .reset_index(
                name="N"
            )
        )

        composition[
            "Species"
        ] = species

        composition[
            "Total_BUSCOs"
        ] = (
            composition
            .groupby(
                "Chromosome"
            )["N"]
            .transform(
                "sum"
            )
        )

        composition[
            "Proportion"
        ] = (
            composition["N"]
            /
            composition["Total_BUSCOs"]
        )

        composition_rows.append(
            composition[
                [
                    "Species",
                    "Chromosome",
                    "Ancestral_element",
                    "N",
                    "Total_BUSCOs",
                    "Proportion",
                ]
            ]
        )

        dominant = assign_chromosomes(
            dat,
            species
        )

        dominant[
            "Assignment_flag"
        ] = "CLEAR"

        dominant.loc[
            dominant[
                "Dominance_margin"
            ] < args.tie_margin,
            "Assignment_flag"
        ] = "NEAR_TIE"

        chromosome_rows.append(
            dominant
        )

        # CLASSIFY EACH BUSCO
        assignment_map = (
            dominant[
                [
                    "Chromosome",
                    "Chromosome_element",
                    "Assignment_flag",
                    "Dominance",
                    "Dominance_margin",
                ]
            ]
        )

        dat = dat.merge(
            assignment_map,
            on="Chromosome",
            how="left"
        )

        dat[
            "Species"
        ] = species

        dat[
            "Observed_element"
        ] = dat[
            "Chromosome_element"
        ]

        dat[
            "Status_movement"
        ] = "Moved"

        dat.loc[
            dat[
                "Ancestral_element"
            ] == dat[
                "Observed_element"
            ],
            "Status_movement"
        ] = "Retained"

        all_gene_rows.append(
            dat[
                [
                    "Species",
                    "BUSCO",
                    "Chromosome",
                    "Ancestral_element",
                    "Observed_element",
                    "Status_movement",
                    "Assignment_flag",
                    "Dominance",
                    "Dominance_margin",
                ]
            ]
        )

        n_chromosomes = (
            dominant[
                "Chromosome"
            ].nunique()
        )

        n_near_ties = (
            dominant[
                "Assignment_flag"
            ] == "NEAR_TIE"
        ).sum()

        print(
            f"Assigned chromosomes: "
            f"{n_chromosomes}"
        )

        print(
            f"Near-tie chromosomes: "
            f"{n_near_ties}"
        )

        print(
            "Chromosomes assigned per Stevens element:"
        )

        element_counts = (
            dominant[
                "Chromosome_element"
            ]
            .value_counts()
        )

        for element, n in (
            element_counts.items()
        ):

            print(
                f"  {element}: {n}"
            )

    if not all_gene_rows:

        raise RuntimeError(
            "No usable gene-movement data were produced."
        )

    genes = pd.concat(
        all_gene_rows,
        ignore_index=True
    )

    chromosomes = pd.concat(
        chromosome_rows,
        ignore_index=True
    )

    composition = pd.concat(
        composition_rows,
        ignore_index=True
    )

    composition.to_csv(
        outdir /
        "chromosome_stevens_composition.tsv",
        sep="\t",
        index=False
    )

    chromosomes.to_csv(
        outdir /
        "chromosome_stevens_assignments.tsv",
        sep="\t",
        index=False
    )

    genes.to_csv(
        outdir /
        "individual_gene_movement.tsv",
        sep="\t",
        index=False
    )

    focal = genes[
        genes[
            "Ancestral_element"
        ].isin(
            FOCAL_ELEMENTS
        )
    ].copy()

    focal.to_csv(
        outdir /
        "individual_gene_movement_XEAGCH.tsv",
        sep="\t",
        index=False
    )

    summary = (
        focal
        .groupby(
            [
                "Species",
                "Ancestral_element"
            ]
        )
        .agg(
            Total_genes=(
                "BUSCO",
                "count"
            ),

            Retained=(
                "Status_movement",
                lambda x:
                    (
                        x == "Retained"
                    ).sum()
            ),

            Moved=(
                "Status_movement",
                lambda x:
                    (
                        x == "Moved"
                    ).sum()
            ),

            Near_tie_genes=(
                "Assignment_flag",
                lambda x:
                    (
                        x == "NEAR_TIE"
                    ).sum()
            ),
        )
        .reset_index()
    )

    summary[
        "Proportion_retained"
    ] = (
        summary[
            "Retained"
        ]
        /
        summary[
            "Total_genes"
        ]
    )

    summary[
        "Proportion_moved"
    ] = (
        summary[
            "Moved"
        ]
        /
        summary[
            "Total_genes"
        ]
    )

    summary[
        "Proportion_near_tie"
    ] = (
        summary[
            "Near_tie_genes"
        ]
        /
        summary[
            "Total_genes"
        ]
    )

    summary.to_csv(
        outdir /
        "gene_retention_summary_XEAGCH.tsv",
        sep="\t",
        index=False
    )

    # MOVEMENT DESTINATIONS

    moved = focal[
        focal[
            "Status_movement"
        ] == "Moved"
    ].copy()

    movement_matrix = (
        moved
        .groupby(
            [
                "Species",
                "Ancestral_element",
                "Observed_element"
            ]
        )
        .size()
        .reset_index(
            name="N_genes"
        )
    )

    movement_matrix.to_csv(
        outdir /
        "gene_movement_destinations_XEAGCH.tsv",
        sep="\t",
        index=False
    )

    # NUMBER OF CHROMOSOMES ASSIGNED TO EACH ELEMENT
    # This helps identify potential fissions.

    element_chromosomes = (
        chromosomes
        .groupby(
            [
                "Species",
                "Chromosome_element"
            ]
        )
        .agg(
            N_chromosomes=(
                "Chromosome",
                "nunique"
            ),

            Chromosomes=(
                "Chromosome",
                lambda x:
                    ",".join(
                        sorted(
                            set(x)
                        )
                    )
            )
        )
        .reset_index()
    )

    element_chromosomes[
        "Possible_fission"
    ] = (
        element_chromosomes[
            "N_chromosomes"
        ] > 1
    )

    element_chromosomes.to_csv(
        outdir /
        "stevens_element_chromosome_counts.tsv",
        sep="\t",
        index=False
    )

    near_ties = (
        chromosomes[
            chromosomes[
                "Assignment_flag"
            ] == "NEAR_TIE"
        ]
        .copy()
    )

    near_ties.to_csv(
        outdir /
        "near_tie_chromosomes.tsv",
        sep="\t",
        index=False
    )

    print()
    print(
        "GENE MOVEMENT ANALYSIS COMPLETE"
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
        "Output files:"
    )

    print(
        "  chromosome_stevens_composition.tsv"
    )

    print(
        "  chromosome_stevens_assignments.tsv"
    )

    print(
        "  individual_gene_movement.tsv"
    )

    print(
        "  individual_gene_movement_XEAGCH.tsv"
    )

    print(
        "  gene_retention_summary_XEAGCH.tsv"
    )

    print(
        "  gene_movement_destinations_XEAGCH.tsv"
    )

    print(
        "  stevens_element_chromosome_counts.tsv"
    )

    print(
        "  near_tie_chromosomes.tsv"
    )


if __name__ == "__main__":
    main()
