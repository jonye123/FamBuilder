#!/bin/bash

# Exit immediately if any command fails.
set -e

# Take in inputs from run_all.sh
PROJECT_NAME="${1:-"myfam"}"
THREADS="${2:-4}"  # This controls how many HMMs run at the same time

# Define the directory containing the .fa files
INPUT_DIR="../${PROJECT_NAME}MSA"     
OUTPUT_DIR="../${PROJECT_NAME}HmmLib"
DB_DIR="../${PROJECT_NAME}HmmDatabase"

# Create output directory and DB directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$DB_DIR"

echo "Building HMMs for ${PROJECT_NAME} using ${THREADS} concurrent jobs..."

# Check for input files
set +e
count=$(ls "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l)
set -e

if [ "$count" -eq 0 ]; then
    echo "Error: No .fasta files found in $INPUT_DIR. Skipping build."
    exit 1
fi

# Define worker function
run_hmmbuild() {
    local file="$1"
    local output_dir="$2"
    
    # Extract HOG ID safely
    local hog_id=$(basename "$file" | grep -oP 'HOG[^_]+')
    local output_file="${output_dir}${hog_id}.hmm"

    # Run HMMER (cpu 1 because 'parallel' handles the concurrency)
    # Use 'if' to check if the command succeeded
    if hmmbuild --cpu 1 -n "$hog_id" "$output_file" "$file" > /dev/null 2>&1; then
        echo "Built HMM for $hog_id"
    else
        echo "Failed to build HMM for $hog_id"
        # We return an error code so parallel knows it failed
        return 1
    fi
}

export -f run_hmmbuild

# get inputs and run the function
find "$INPUT_DIR" -name "*.fasta" | parallel --gnu -j "$THREADS" run_hmmbuild {} "$OUTPUT_DIR"

echo "All HMMs have been built."