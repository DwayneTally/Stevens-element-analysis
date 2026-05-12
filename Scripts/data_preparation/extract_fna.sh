#!/bin/bash

#Directory to store extracted .fna files
mkdir -p ../genomic_fna
#Loop over all genome ZIP files
for zip_file in *_genome.zip; do
    #Extract the species and assembly name
    base_name="${zip_file%_genome.zip}"
    output_file="../genomic_fna/${base_name}_genomic.fna"
    fna_path=$(unzip -l "$zip_file" | awk '{print $4}' | grep -m 1 '\.fna$')
    if [ -n "$fna_path" ]; then
        #Extract the .fna file from the ZIP archive
        unzip -p "$zip_file" "$fna_path" > "$output_file"
        echo "Extracted $output_file"
    else
        echo "Warning: No .fna file found in $zip_file"
    fi
done

echo "All .fna files have been extracted to the genomic_fna directory."

