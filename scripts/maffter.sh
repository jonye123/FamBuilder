#!/bin/bash
set -e

# maffter.sh - Run MAFFT alignment on HOG FASTA files
# Skips HOGs whose FASTA hasn't changed since the previous iteration (MSA caching)
#
# Usage: bash maffter.sh --input_dir <path> --output_dir <path> --jobs <n> [--prev_input_dir <path>] [--prev_msa_dir <path>]

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

input_dir=""
output_dir=""
jobs=8
prev_input_dir=""
prev_msa_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input_dir)    input_dir="$2";     shift 2 ;;
        --output_dir)   output_dir="$2";    shift 2 ;;
        --jobs)         jobs="$2";          shift 2 ;;
        --prev_input_dir) prev_input_dir="$2"; shift 2 ;;
        --prev_msa_dir)   prev_msa_dir="$2";   shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "$input_dir" || -z "$output_dir" ]]; then
    echo "Error: --input_dir and --output_dir are required."
    exit 1
fi

if [[ ! -d "$input_dir" ]]; then
    echo "Error: input_dir '$input_dir' does not exist."
    exit 1
fi

mkdir -p "$output_dir"

echo "Running MAFFT alignments..."
echo "  Input dir  : $input_dir"
echo "  Output dir : $output_dir"
echo "  Jobs       : $jobs"
echo "  Prev input : ${prev_input_dir:-none}"
echo "  Prev MSA   : ${prev_msa_dir:-none}"

# ---------------------------------------------------------------------------
# MSA caching: compare checksums to detect unchanged HOGs
# ---------------------------------------------------------------------------

# Build a list of FASTA files that actually need realigning
to_align=()
copied=0
to_run=0

for fasta_file in "$input_dir"/*.fasta; do
    [[ -f "$fasta_file" ]] || continue

    filename=$(basename -- "$fasta_file" .fasta)
    output_file="${output_dir}/${filename}_mafft_output.fasta"
    cache_hit=false

    # Check if we have a previous iteration to compare against
    if [[ -n "$prev_input_dir" && -n "$prev_msa_dir" ]]; then
        prev_fasta="${prev_input_dir}/${filename}.fasta"
        prev_msa="${prev_msa_dir}/${filename}_mafft_output.fasta"

        if [[ -f "$prev_fasta" && -f "$prev_msa" ]]; then
            # Compare checksums
            curr_checksum=$(md5sum "$fasta_file" | awk '{print $1}')
            prev_checksum=$(md5sum "$prev_fasta" | awk '{print $1}')

            if [[ "$curr_checksum" == "$prev_checksum" ]]; then
                # FASTA unchanged — reuse previous MSA
                cp "$prev_msa" "$output_file"
                cache_hit=true
                ((copied++)) || true
            fi
        fi
    fi

    if [[ "$cache_hit" == false ]]; then
        to_align+=("$fasta_file")
        ((to_run++)) || true
    fi
done

echo "  Cache hits (MSA reused) : $copied"
echo "  To align                : $to_run"

# ---------------------------------------------------------------------------
# MAFFT alignment (parallelized)
# ---------------------------------------------------------------------------

mafft_file() {
    local file="$1"
    local output_dir="$2"
    local filename=$(basename -- "$file" .fasta)
    local output_file="${output_dir}/${filename}_mafft_output.fasta"

    mafft --localpair --maxiterate 1000 --anysymbol "$file" > "$output_file"
    if [[ $? -ne 0 ]]; then
        echo "Error: MAFFT failed for '$file'." >&2
        exit 1
    fi
}

export -f mafft_file

if [[ ${#to_align[@]} -gt 0 ]]; then
    printf '%s\n' "${to_align[@]}" | parallel --gnu -j "$jobs" mafft_file {} "$output_dir"
else
    echo "All MSAs were cached — no MAFFT runs needed."
fi

echo "MAFFT alignments complete."