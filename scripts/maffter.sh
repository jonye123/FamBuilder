#!/bin/bash
# Usage: ./maffter.sh [project_name] [threads]

# Exit immediately if any command fails.
set -e

# Take in inputs from run_all.sh
PROJECT_NAME=${1:-"myfam"}
THREADS=${2:-4}

# Define the directory containing the .fa files
# These are relative to where the script is run
INPUT_DIR="hogDirectoryFiltered"
OUTPUT_DIR="hogMsaDir/"
# OUTPUT_DIR="../${PROJECT_NAME}MSA/"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Define the function for GNU Parallel
mafft_file() {
    local file="$1"
    local output_dir="$2"
    
    local filename=$(basename -- "$file" .fasta)
    local output_file="${output_dir}${filename}_mafft_output.fasta"
    
    # Run MAFFT (Quiet mode to reduce log spam)
    mafft --quiet --localpair --maxiterate 1000 --anysymbol "$file" > "$output_file"
}

export -f mafft_file

echo "Running MAFFT alignment with $THREADS threads..."

# Run Parallel
# We pass output_dir as the second argument to the function
find "$INPUT_DIR" -type f -name '*.fasta' | parallel --gnu -j "$THREADS" mafft_file {} "$OUTPUT_DIR"

echo "MAFFT alignment complete."