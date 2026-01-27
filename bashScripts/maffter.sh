#!/bin/bash

# Define the directory containing the .fa files
input_dir="../hogDirectoryFiltered"
output_dir="../MSA/"

mkdir -p $output_dir

mafft_file() {
    local file="$1"
    local filename=$(basename -- "$file" .fasta)
    local output_file="${output_dir}${filename}_mafft_output.fasta"
    mafft --localpair --maxiterate 1000 --anysymbol "$file" > "$output_file"
}

export -f mafft_file
export output_dir

find "$input_dir" -type f -name '*.fasta' | parallel --gnu -j 16 mafft_file {} {/}