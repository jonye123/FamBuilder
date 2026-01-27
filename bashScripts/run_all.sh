#!/bin/bash

log_dir="bashLogs"
mkdir "$log_dir"

eval "$(micromamba shell hook bash)"
micromamba activate orthofinder

bash maffter.sh > "$log_dir/maffter.log" 2> "$log_dir/maffter.err"

bash hmmbuilder.sh > "$log_dir/grassHmm.log" 2> "$log_dir/grassHmm.err"

cat ../grassfamHmmLib/* > ../grassHmmDatabase/grassfam 
hmmpress ../grassHmmDatabase/grassfam

bash confusionMatrixMaker.sh > "$log_dir/confusion.log" 2> "$log_dir/confusion.err"

wait
echo "All processes finished."