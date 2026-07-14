#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# confusionMatrixMaker.sh
# Run hmmscan on each HOG's sequences against the full HMM library and
# produce a confusion matrix for evaluation.
# Part of the FamBuilder pipeline.
#
# Usage:
#   bash confusionMatrixMaker.sh \
#       --hog_dir <hogDirectoryFiltered> \
#       --hmm_lib <path/to/pressedHmmLib> \
#       --output_dir <database_dir> \
#       --jobs <n>
# ---------------------------------------------------------------------------

# --- Argument parsing ---
hog_dir=""
hmm_lib=""
output_dir=""
jobs=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hog_dir)    hog_dir="$2";    shift 2 ;;
        --hmm_lib)    hmm_lib="$2";    shift 2 ;;
        --output_dir) output_dir="$2"; shift 2 ;;
        --jobs)       jobs="$2";       shift 2 ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash confusionMatrixMaker.sh --hog_dir <dir> --hmm_lib <lib> --output_dir <dir> --jobs <n>"
            exit 1
            ;;
    esac
done

# --- Validate arguments ---
if [[ -z "$hog_dir" || -z "$hmm_lib" || -z "$output_dir" ]]; then
    echo "Error: --hog_dir, --hmm_lib, and --output_dir are all required."
    echo "Usage: bash confusionMatrixMaker.sh --hog_dir <dir> --hmm_lib <lib> --output_dir <dir> --jobs <n>"
    exit 1
fi

if [[ ! -d "$hog_dir" ]]; then
    echo "Error: HOG directory '$hog_dir' does not exist."
    exit 1
fi

if [[ ! -f "$hmm_lib" ]]; then
    echo "Error: HMM library '$hmm_lib' does not exist."
    exit 1
fi

# --- Check dependencies ---
if ! command -v hmmscan &> /dev/null; then
    echo "Error: hmmscan not found in PATH."
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found in PATH."
    exit 1
fi

# --- Set up output directories ---
scan_dir="${output_dir}/hmmscan_results"
cat_dir="${output_dir}/catThis"
confusion_path="${output_dir}/confusionMatrix"

mkdir -p "$scan_dir"
mkdir -p "$cat_dir"

# Clear output files to avoid double-append on rerun
> "$confusion_path"
# Write header
echo -e "Gene\tHog\tHmmHogHit\tE-value\tScore" > "$confusion_path"

# Clear cat directory to avoid stale results
rm -f "$cat_dir"/*

echo "Running hmmscan on HOGs in '$hog_dir' against '$hmm_lib'..."
echo "Output going to '$output_dir'..."

# --- hmmscan function ---
process_file() {
    local file="$1"
    local hmm_lib="$2"
    local scan_dir="$3"
    local cat_dir="$4"

    # Extract HOG ID robustly:
    # Handles N0.HOG001.fasta, HOG001.fasta, and suffixed N0.HOG001-1.fasta
    local basename
    basename=$(basename "$file" .fasta)
    # Remove N0. prefix if present
    basename="${basename#N0.}"
    local hog_id="$basename"

    if [[ -z "$hog_id" ]]; then
        echo "Warning: Could not extract HOG ID from '$file'. Skipping."
        return 1
    fi

    local scan_output="${scan_dir}/${hog_id}"

    # Run hmmscan
    if ! hmmscan --tblout "$scan_output" "$hmm_lib" "$file" &> /dev/null; then
        echo "Error: hmmscan failed for $hog_id."
        return 1
    fi

    # Parse results: for each gene take only the top hit (first non-comment line per gene)
    awk -v hog_id="$hog_id" \
        '!/^#/ && !seen[$3]++ {print $3"\t"hog_id"\t"$1"\t"$5"\t"$6}' \
        "$scan_output" > "${cat_dir}/${hog_id}"
}

export -f process_file

# --- Run in parallel ---
find "$hog_dir" -type f -name "*.fasta" | \
    parallel --gnu \
             -j "$jobs" \
             --halt soon,fail=1 \
             --joblog "${output_dir}/hmmscan.joblog" \
             process_file {} "$hmm_lib" "$scan_dir" "$cat_dir"

# --- Concatenate results into confusion matrix ---
# Check that we have results to concatenate
if [[ -z "$(ls -A "$cat_dir")" ]]; then
    echo "Error: No hmmscan results found in '$cat_dir'. Something went wrong."
    exit 1
fi

find "$cat_dir" -type f -print0 | xargs -0 cat >> "$confusion_path"

echo "Confusion matrix written to '$confusion_path'."
echo "Done."