#!/bin/bash

# GENERATE GENE PRESENCE/ABSENCE MATRIX
# This script reads filtered BLAST hits for each sample,
# checks which target genes are present, and generates
# a tab-separated presence/absence matrix.

# Directories and files
INPUT_DIR="./blast_best_hits"        # Filtered hits from previous step
TARGET_GENES="./target_genes.txt"    # List of target genes
OUTPUT_MATRIX="./gene_presence_absence_matrix.tsv"  # Output matrix

# Read target genes into array
mapfile -t genes < "$TARGET_GENES"

# Prepare header row
{
    echo -ne "Sample"
    for g in "${genes[@]}"; do
        echo -ne "\t$g"
    done
    echo
} > "$OUTPUT_MATRIX"

# Loop over each sample file
for f in "$INPUT_DIR"/*.tsv; do
    sample=$(basename "$f" .tsv)
    
    # Initialize gene presence hash
    declare -A present
    for g in "${genes[@]}"; do
        present["$g"]=0
    done

    # Mark present genes from BLAST hits
    while IFS=$'\t' read -r _ hit_gene _; do
        # Extract gene name from column 2 (before first "_")
        gname=$(echo "$hit_gene" | cut -d'_' -f1)
        if [[ -n "${present[$gname]}" ]]; then
            present[$gname]=1
        fi
    done < "$f"

    # Print sample row
    echo -ne "$sample"
    for g in "${genes[@]}"; do
        echo -ne "\t${present[$g]}"
    done
    echo
done >> "$OUTPUT_MATRIX"

echo "Presence/absence matrix saved to $OUTPUT_MATRIX"