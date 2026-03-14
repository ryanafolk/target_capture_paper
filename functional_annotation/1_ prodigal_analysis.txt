```bash
#!/bin/bash

# Run Prodigal gene prediction on all metagenomic contigs
# This script predicts genes from assembled contigs using Prodigal in metagenome mode

# Directory containing contig FASTA files
CONTIGS_DIR="./contigs"

# Directory where Prodigal output files will be saved
OUTPUT_DIR="./prodigal_outputs"

# Create the output directory if it does not already exist
mkdir -p "$OUTPUT_DIR"

# Loop through every FASTA file in the contigs directory
for contig_file in "$CONTIGS_DIR"/*.fasta; do

    # Extract the filename without the directory path and remove .fasta extension
    sample=$(basename "$contig_file" .fasta)

    # Remove the "_scaffolds" suffix to keep a clean sample name
    sample=${sample%_scaffolds}

    # Print sample name to track progress
    echo "Processing sample: $sample"

    # Create a separate output folder for each sample
    sample_dir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_dir"

    # Run Prodigal in metagenome mode
    # -i : input contig file
    # -o : gene predictions in GenBank format
    # -a : predicted protein sequences (FAA)
    # -d : nucleotide sequences of predicted genes (FNA)
    # -p meta : metagenome mode for fragmented assemblies
    prodigal -i "$contig_file" \
             -o "$sample_dir/${sample}_genes.gbk" \
             -a "$sample_dir/${sample}_proteins.faa" \
             -d "$sample_dir/${sample}_genes.fna" \
             -p meta

done

# Message indicating completion of gene prediction for all samples
echo "Prodigal analysis complete for all samples."
```
