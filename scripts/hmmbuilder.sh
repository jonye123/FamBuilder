#!/bin/bash

# Define the directory containing the .fa files
input_dir="../MSA/"
output_dir="../grassfamHmmLib/"
for_later="../grassHmmDatabase/"

mkdir -p "$output_dir"

mkdir -p "$for_later"

# Array to store background process IDs and exit statuses
declare -A hmmbuild_exit_status

# Iterate over each .fa file in the directory
for file in "$input_dir"*.fasta; do
    if [ -f "$file" ] && [ -r "$file" ]; then
        # Get the base filename without the directory path and extension
        hog_id=$(basename "$file" | grep -oP 'HOG[^_]+')
        # Define the output file name
        output_file="${output_dir}${hog_id}.hmm"
        # Generate the hmm command and run it in background
        hmmbuild -n "$hog_id" --cpu "$(nproc)" "$output_file" "$file" &
        # Store the PID of the background process
        pid=$!
        echo "pid: ${pid}, hog_id: ${hog_id}"
        # Store the PID and its exit status
        hmmbuild_exit_status["$pid"]=$?
    else
        echo "Skipped file '$file'. It does not exist or is not readable."
    fi
done

# Wait for all background processes to finish
wait

# Check the exit status of each background process
for pid in "${!hmmbuild_exit_status[@]}"; do
    if [ "${hmmbuild_exit_status[$pid]}" -eq 0 ]; then
        echo "Built HMM for ${input_dir}${pid}"
    else
        echo "Failed to build HMM for ${input_dir}${pid}"
    fi
done

echo "All HMMs have been built."