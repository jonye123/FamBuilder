#!/bin/bash
# Usage: ./confusionMatrixMaker.sh [project_name] [threads]
set -e

# Capture Inputs
PROJECT_NAME="${1:-"myfam"}"
THREADS="${2:-4}"

# Define Generic Directories
# Input: Where the filtered sequences are
SOURCE_DIRECTORY="../hogDirectoryFiltered"
# Base Database Path (Matches main.sh expectation: ProjectHmmDatabase/Projectfam)
BASE_DB_PATH="../${PROJECT_NAME}HmmDatabase/${PROJECT_NAME}fam"

# Specific paths inside the DB folder
DESTINATION_DIRECTORY="${BASE_DB_PATH}/test"
CONFUSION_DIRECTORY="${BASE_DB_PATH}/catThis"
CONFUSION_PATH="${BASE_DB_PATH}/confusionMatrix"

# Create new directories
mkdir -p "$DESTINATION_DIRECTORY"
mkdir -p "$CONFUSION_DIRECTORY"

# The HMM Database pressed in run_all.sh after building the HMM
HMM_DB_PATH="../${PROJECT_NAME}HmmDatabase/${PROJECT_NAME}fam1"

# Initialize the Confusion Matrix File
echo -e "Gene\tHog\tHmmHogHit\tE-value\tScore" > "$CONFUSION_PATH"

echo "Calculating Confusion Matrix for $PROJECT_NAME..."

# Define the Worker Function
process_file() {
    local file="$1"
    local hog_id=$(basename "$file" | grep -oP 'HOG[^.]+') # Extract HOG ID
    local output_file="${DESTINATION_DIRECTORY}/${hog_id}"
    
    # Run hmmscan
    # We scan the gene file ($file) against our HMM database ($HMM_DB_PATH)
    hmmscan --tblout "$output_file" "$HMM_DB_PATH" "$file" &> /dev/null
    
    # Parse the output
    # We take the best hit for each query and write to a temp file
    awk -v hog_id="$hog_id" '!/^#/ && !seen[$3]++ {print $3"\t"hog_id"\t"$1"\t"$5"\t"$6}' "$output_file" > "$CONFUSION_DIRECTORY/$hog_id"
}

export -f process_file
export DESTINATION_DIRECTORY
export CONFUSION_DIRECTORY
export HMM_DB_PATH

find "$SOURCE_DIRECTORY" -type f -name '*.fasta' | parallel --gnu -j $THREADS process_file {} {/}

cat "$confusion_directory"/* >> "$confusion_path"