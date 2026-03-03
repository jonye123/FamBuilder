#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# hmmbuilder.sh
# Build HMMs from multiple sequence alignments using hmmbuild.
# Part of the FamBuilder pipeline.
#
# Usage:
#   bash hmmbuilder.sh --input_dir <MSA_dir> --output_dir <HmmLib_dir> --jobs <n>
# ---------------------------------------------------------------------------

# --- Argument parsing ---
input_dir=""
output_dir=""
jobs=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input_dir)  input_dir="$2";  shift 2 ;;
        --output_dir) output_dir="$2"; shift 2 ;;
        --jobs)       jobs="$2";       shift 2 ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash hmmbuilder.sh --input_dir <MSA_dir> --output_dir <HmmLib_dir> --jobs <n>"
            exit 1
            ;;
    esac
done

# --- Validate arguments ---
if [[ -z "$input_dir" || -z "$output_dir" ]]; then
    echo "Error: --input_dir and --output_dir are required."
    echo "Usage: bash hmmbuilder.sh --input_dir <MSA_dir> --output_dir <HmmLib_dir> --jobs <n>"
    exit 1
fi

if [[ ! -d "$input_dir" ]]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

mkdir -p "$output_dir"

# --- Check dependencies ---
if ! command -v hmmbuild &> /dev/null; then
    echo "Error: hmmbuild not found in PATH."
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found in PATH."
    exit 1
fi

echo "Building HMMs from '$input_dir' -> '$output_dir' using $jobs jobs..."

# --- HMM build function ---
build_hmm() {
    local file="$1"
    local output_dir="$2"

    # Extract HOG ID robustly:
    # Handles N0.HOG001_mafft_output.fasta, HOG001_mafft_output.fasta,
    # and suffixed files like N0.HOG001-1_mafft_output.fasta
    local basename
    basename=$(basename "$file" .fasta)
    # Remove _mafft_output suffix if present
    basename="${basename%_mafft_output}"
    # Remove N0. prefix if present
    basename="${basename#N0.}"
    local hog_id="$basename"

    if [[ -z "$hog_id" ]]; then
        echo "Warning: Could not extract HOG ID from '$file'. Skipping."
        return 1
    fi

    local output_file="${output_dir}/${hog_id}.hmm"

    echo "Building HMM for $hog_id..."
    if hmmbuild -n "$hog_id" --cpu 1 "$output_file" "$file" &> "${output_dir}/${hog_id}.log"; then
        echo "Done: $hog_id"
        return 0
    else
        echo "Error: hmmbuild failed for $hog_id. Check ${output_dir}/${hog_id}.log"
        return 1
    fi
}

export -f build_hmm
export output_dir

# --- Run in parallel using GNU parallel ---
# --halt soon,fail=1 stops all jobs if any fail
# --joblog tracks which jobs succeeded/failed
find "$input_dir" -type f -name "*.fasta" | \
    parallel --gnu \
             -j "$jobs" \
             --halt soon,fail=1 \
             --joblog "${output_dir}/hmmbuild.joblog" \
             build_hmm {} "$output_dir"

echo "All HMMs built successfully. Output in '$output_dir'."