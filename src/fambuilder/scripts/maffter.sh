#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# maffter.sh
# Run MAFFT alignments on HOG FASTA files.
# Implements MSA caching — if a HOG's FASTA hasn't changed since the previous
# iteration, the existing MSA is copied instead of rerunning MAFFT.
#
# Part of the FamBuilder pipeline.
#
# Usage:
#   bash maffter.sh \
#       --input_dir <hogDirectoryFiltered> \
#       --output_dir <MSA_dir> \
#       --jobs <n> \
#       [--prev_msa_dir <previous_MSA_dir>] \
#       [--prev_checksum_file <previous_checksum_file>]
# ---------------------------------------------------------------------------

# --- Argument parsing ---
input_dir=""
output_dir=""
jobs=8
prev_msa_dir=""
prev_checksum_file=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input_dir)           input_dir="$2";           shift 2 ;;
        --output_dir)          output_dir="$2";          shift 2 ;;
        --jobs)                jobs="$2";                shift 2 ;;
        --prev_msa_dir)        prev_msa_dir="$2";        shift 2 ;;
        --prev_checksum_file)  prev_checksum_file="$2";  shift 2 ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash maffter.sh --input_dir <dir> --output_dir <dir> --jobs <n> [--prev_msa_dir <dir>] [--prev_checksum_file <file>]"
            exit 1
            ;;
    esac
done

# --- Validate required arguments ---
if [[ -z "$input_dir" || -z "$output_dir" ]]; then
    echo "Error: --input_dir and --output_dir are required."
    exit 1
fi

if [[ ! -d "$input_dir" ]]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# --- Check dependencies ---
if ! command -v mafft &> /dev/null; then
    echo "Error: mafft not found in PATH."
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found in PATH."
    exit 1
fi

if ! command -v md5sum &> /dev/null; then
    echo "Error: md5sum not found in PATH."
    exit 1
fi

mkdir -p "$output_dir"

# --- Helper: look up previous checksum for a HOG ID using grep ---
get_prev_checksum() {
    local hog_id="$1"
    if [[ -n "$prev_checksum_file" && -f "$prev_checksum_file" ]]; then
        grep -m1 "^${hog_id}	" "$prev_checksum_file" | cut -f2
    fi
}

# --- Current checksum file for this iteration ---
checksum_file="${output_dir}/msa_checksums.tsv"
> "$checksum_file"

# --- Counters ---
count_cached=0
count_rerun=0

echo "--- FamBuilder maffter.sh ---"
echo "Input dir : $input_dir"
echo "Output dir: $output_dir"
echo "Jobs      : $jobs"

# --- Separate HOGs into cached and needs-rerun ---
cached_list=()
rerun_list=()

for file in "$input_dir"/*.fasta; do
    [[ -f "$file" ]] || continue

    # Extract HOG ID — strip N0. prefix
    basename=$(basename "$file" .fasta)
    hog_id="${basename#N0.}"

    # Compute checksum of current FASTA
    current_checksum=$(md5sum "$file" | awk '{print $1}')

    # Write to current checksum file
    echo -e "${hog_id}\t${current_checksum}" >> "$checksum_file"

    # Determine expected MSA output filename
    prev_msa="${prev_msa_dir}/${basename}_mafft_output.fasta"

    # Check if we can use the cached MSA
    prev_checksum=$(get_prev_checksum "$hog_id")
    if [[ -n "$prev_checksum" ]] && \
       [[ "$prev_checksum" == "$current_checksum" ]] && \
       [[ -n "$prev_msa_dir" ]] && \
       [[ -f "$prev_msa" ]]; then
        cached_list+=("$file")
    else
        rerun_list+=("$file")
    fi
done

echo "HOGs to rerun MAFFT : ${#rerun_list[@]}"
echo "HOGs to copy (cache): ${#cached_list[@]}"

# --- Copy cached MSAs ---
if [[ ${#cached_list[@]} -gt 0 ]]; then
    echo "Copying cached MSAs..."
    for file in "${cached_list[@]}"; do
        basename=$(basename "$file" .fasta)
        prev_msa="${prev_msa_dir}/${basename}_mafft_output.fasta"
        cp "$prev_msa" "${output_dir}/${basename}_mafft_output.fasta"
        count_cached=$((count_cached + 1))
    done
    echo "Copied $count_cached cached MSAs."
fi

# --- Run MAFFT on changed/new HOGs ---
if [[ ${#rerun_list[@]} -gt 0 ]]; then
    echo "Running MAFFT on ${#rerun_list[@]} HOGs..."

    mafft_file() {
        local file="$1"
        local output_dir="$2"
        local filename
        filename=$(basename "$file" .fasta)
        local output_file="${output_dir}/${filename}_mafft_output.fasta"

        if mafft --auto --anysymbol "$file" > "$output_file" 2>/dev/null; then
            echo "Done: $filename"
        else
            echo "Error: MAFFT failed for $filename"
            return 1
        fi
    }

    export -f mafft_file

    printf '%s\n' "${rerun_list[@]}" | \
        parallel --gnu \
                 -j "$jobs" \
                 --halt soon,fail=1 \
                 --joblog "${output_dir}/mafft.joblog" \
                 mafft_file {} "$output_dir"

    count_rerun=${#rerun_list[@]}
    echo "MAFFT complete for $count_rerun HOGs."
fi

echo ""
echo "MSA checksums written to '$checksum_file'."
echo "--- maffter.sh finished: $count_cached cached, $count_rerun rerun ---"