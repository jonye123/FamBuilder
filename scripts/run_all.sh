#!/bin/bash

# Usage: run_all.sh [project_name] [threads]

# Exit immediately if any command fails.
set -e

# Take in inputs
PROJECT_NAME=${1:-"myfam"}
THREADS=${2:-4}

# Define generic paths
LOG_DIR="bashLogs"
HMM_LIB_DIR="${PROJECT_NAME}HmmLib"
DB_DIR="${PROJECT_NAME}HmmDatabase"

# Create log directory
mkdir -p "$LOG_DIR"

# Assume the environment is already activated by main.sh or the user.
# eval "$(micromamba shell hook bash)"
# micromamba activate orthofinder

echo "Starting processes for project: $PROJECT_NAME..."

# Run MAFFT (alignment) by calling maft script
bash maffter.sh ${PROJECT_NAME} ${THREADS} > "$LOG_DIR/maffter.log" 2> "$LOG_DIR/maffter.err"

# Run HMM Builder
bash hmmbuilder.sh ${PROJECT_NAME} ${THREADS} > "$LOG_DIR/HmmBuilder.log" 2> "$LOG_DIR/HmmBuilder.err"

# Press the HMMs (Prepare the database)
# Ensure the DB directory exists
mkdir -p "$DB_DIR"

# Check if the HMM library actually has files before trying to cat them
if ls "$HMM_LIB_DIR"/*.hmm 1> /dev/null 2>&1; then
    cat "$HMM_LIB_DIR"/*.hmm > "$DB_DIR/${PROJECT_NAME}fam1"
    hmmpress -f "$DB_DIR/${PROJECT_NAME}fam1"
else
    echo "Error: No .hmm files found in $HMM_LIB_DIR"
    exit 1
fi

# Run Confusion Matrix
bash confusionMatrixMaker.sh ${PROJECT_NAME} ${THREADS} > "$LOG_DIR/confusion.log" 2> "$LOG_DIR/confusion.err"

wait
echo "All processes finished."