#!/usr/bin/env python3

import subprocess
from pathlib import Path

# === PATH SETUP ===
gffread_bin = Path("/N/project/Bracewell_beetle/SOFTWARE/gffread/gffread/gffread")
gff_dir = Path("Gabby_genomes")
genome_dir = Path("Gabby_genomes")  # location of unmasked .fna files
pep_output_dir = Path("Gabby_peptide")
cds_output_dir = Path("Gabby_cds")

# === CREATE OUTPUT DIRS ===
pep_output_dir.mkdir(exist_ok=True)
cds_output_dir.mkdir(exist_ok=True)

# === LOOP OVER GFF3 FILES ===
for gff_file in gff_dir.glob("*.gff"):
    base = gff_file.stem

    # Point to the unmasked genome file
    genome_path = genome_dir / f"{base}_genomic.fna"
    if not genome_path.exists():
        print(f"⚠️  Skipping {base} — genome not found at {genome_path}")
        continue

    # Output file paths
    cds_out = cds_output_dir / f"{base}_cds.fa"
    pep_out = pep_output_dir / f"{base}.fa"

    # Run gffread
    cmd = [
        str(gffread_bin),
        str(gff_file),
        "-g", str(genome_path),
        "-x", str(cds_out),
        "-y", str(pep_out)
    ]

    print(f"🔁 Running gffread for {base}")
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Done: {base}")
    except subprocess.CalledProcessError:
        print(f"❌ Failed: {base}")

