#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# splitAndScore.sh
# Split a (large, per-species) peptide FASTA into chunks and hmmscan each
# chunk against an HMM library in parallel, producing per-species confusion
# rows. Better suited than confusionMatrixMaker.sh for scoring a small number
# of large FASTAs (e.g. whole proteomes) rather than many small HOG FASTAs.
#
# If hmmscan errors or segfaults partway through a chunk, the chunk is
# trimmed to just after the last successfully scored gene and retried,
# so a single bad sequence doesn't lose progress on the rest of the file.
#
# Part of the FamBuilder pipeline.
#
# Requires: hmmscan, seqkit, GNU parallel (all in PATH)
#
# Usage:
#   bash splitAndScore.sh \
#       --source_fasta <path/to/peptide.fasta> \
#       --hmm_lib <path/to/pressed/hmm/lib> \
#       --jobs <n> \
#       [--split_dir <dir>] \
#       [--output_dir <dir>] \
#       [--confusion_dir <dir>]
# ---------------------------------------------------------------------------

# --- Argument parsing ---
source_fasta=""
hmm_lib=""
jobs=8
split_dir="./split"
output_dir="./hmmerOutput"
confusion_dir="./catThis"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source_fasta) source_fasta="$2"; shift 2 ;;
        --hmm_lib)      hmm_lib="$2";      shift 2 ;;
        --jobs)         jobs="$2";         shift 2 ;;
        --split_dir)    split_dir="$2";    shift 2 ;;
        --output_dir)   output_dir="$2";   shift 2 ;;
        --confusion_dir) confusion_dir="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash splitAndScore.sh --source_fasta <fasta> --hmm_lib <lib> --jobs <n> [--split_dir <dir>] [--output_dir <dir>] [--confusion_dir <dir>]"
            exit 1
            ;;
    esac
done

# --- Validate arguments ---
if [[ -z "$source_fasta" || -z "$hmm_lib" ]]; then
    echo "Error: --source_fasta and --hmm_lib are required."
    echo "Usage: bash splitAndScore.sh --source_fasta <fasta> --hmm_lib <lib> --jobs <n> [--split_dir <dir>] [--output_dir <dir>] [--confusion_dir <dir>]"
    exit 1
fi

if [[ ! -f "$source_fasta" ]]; then
    echo "Error: Source FASTA '$source_fasta' does not exist."
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

if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit not found in PATH."
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found in PATH."
    exit 1
fi

start_time=$(date +%s)

mkdir -p "$output_dir"
mkdir -p "$confusion_dir"
mkdir -p "$split_dir"
mkdir -p ./hmmerLogs

# --- Per-chunk scoring function ---
process_file() {
    local file="$1"
    local species_name
    species_name=$(basename "$file")
    local output_file="${output_dir}/${species_name}"
    local letter
    letter=$(echo "$species_name" | grep -o '^[a-zA-Z]*')

    : > "$confusion_dir/$species_name" # empty the confusion file if it exists
    : > ./hmmerLogs/"$species_name".log # empty the log file if it exists
    : > "$output_file"

    echo "Processing: $species_name"

    while [[ -s "$file" ]]; do
        tmp_tblout=$(mktemp /tmp/tmp_tblout.${species_name}.XXXXXX)
        hmmscan --cpu 1 --tblout "$tmp_tblout" "$hmm_lib" "$file" >> ./hmmerLogs/"$species_name".log 2>&1

        local status=$?

        awk -v species_name="$species_name" -v letter="$letter" '!/^#/ && !seen[$3]++ {print $3"\t"letter"\t"$1"\t"$5"\t"$6}' "$tmp_tblout" >> "$confusion_dir/$species_name"
        cat "$tmp_tblout" >> "$output_file"
        rm -f "$tmp_tblout"

        if [[ $status -eq 0 ]]; then
            echo "hmmscan completed successfully for $species_name"
            break
        else
            echo "Segfault or error ($status) detected for $species_name, trimming fasta and retrying..."

            # Get last successfully parsed gene
            last_gene=$(tail -n 1 "$confusion_dir/$species_name" | cut -f1)

            # Subset the fasta to continue from the next next gene (since the gene after the last parsed one is likely the issue)
            if [[ -n "$last_gene" ]]; then
                echo "Last successfully parsed gene: $last_gene"
                # Subset fasta to start AFTER last_gene (skip the immediate next gene too)
                awk -v last="$last_gene" 'BEGIN { found=0; skip=0 } /^>/ { if ($0 ~ ">"last"$") { found=1; skip=2 } else if (found && skip > 0) { skip-- } else if (found && skip == 0) { printing=1 } } { if (printing) print }' "$file" > "${file}.trimmed"

                mv "${file}.trimmed" "$file"
            else
                echo "No progress possible on $species_name, aborting."
                break
            fi
        fi
    done
}

export -f process_file
export output_dir
export hmm_lib
export confusion_dir

seqkit split -p "$jobs" -O "$split_dir" "$source_fasta"

find "$split_dir" -type f -name '*' | parallel --gnu -j "$jobs" process_file {}

end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Confusion rows written to '$confusion_dir'."
echo "Total execution time: $execution_time seconds"
