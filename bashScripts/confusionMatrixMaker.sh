#!/bin/bash

source_directory="../hogDirectoryFiltered"
destination_directory="../grassHmmDatabase/grassfam1/test1"
grassfam_path="../grassHmmDatabase/grassfam"
confusion_path="../grassHmmDatabase/grassfam1/confusionMatrix"
confusion_directory="../grassHmmDatabase/grassfam1/catThis"

mkdir -p "$destination_directory"
mkdir -p "$confusion_directory"

echo -e "Gene\tHog\tHmmHogHit\tE-value\tScore" > "$confusion_path"

process_file() {
    local file="$1"
    local hog_id=$(basename "$file" | grep -oP 'HOG[^.]+')
    local output_file="${destination_directory}/${hog_id}"
    echo "$hog_id"
    hmmscan --tblout "$output_file" "$grassfam_path" "$file" &> /dev/null
    awk -v hog_id="$hog_id" '!/^#/ && !seen[$3]++ {print $3"\t"hog_id"\t"$1"\t"$5"\t"$6}' "$output_file" > "$confusion_directory/$hog_id"
}

export -f process_file
export destination_directory
export grassfam_path
export confusion_directory

find "$source_directory" -type f -name '*.fasta' | parallel --gnu -j 16 process_file {} {/}

cat "$confusion_directory"/* >> "$confusion_path"