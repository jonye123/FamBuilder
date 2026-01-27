#!/bin/bash

# GrassFam Iterative HMM-Based Clustering Algorithm
# Input: N0.tsv (gene families from OrthoFinder)

input_file="N0.tsv"
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file $input_file not found."
    exit 1
fi

eval "$(micromamba shell hook bash)"
micromamba activate orthofinder

# Initialize iteration counter
iteration=1
prev_error=""

while true; do
    dir="HmmRun${iteration}"

    echo "Starting iteration $iteration in directory $dir"

    # Create the directory
    mkdir -p "$dir"/bashScripts

    # Copy necessary files
    cp bashScripts/* "$dir"/bashScripts
    cp confusionMatrxPlots.ipynb "$dir"/
    cp hogFilter.ipynb "$dir"/
    cp orthofinder_to_fasta.py "$dir"/
    cp "$input_file" "$dir"

    # Move into the directory
    cd "$dir" || exit 1

    # Execute hogFilter notebook
    jupyter nbconvert --to notebook --execute --inplace hogFilter.ipynb
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to execute hogFilter.ipynb in iteration $iteration."
        exit 1
    fi

    # Move into bashScripts directory

    cd bashScripts || exit 1

    # Run the bash script
    bash run_all.sh &> all.log
    if [[ $? -ne 0 ]]; then
        echo "Error: run_all.sh failed in iteration $iteration."
        exit 1
    fi

    # Move back to the parent directory
    cd ..

    # Execute confusionMatrxPlots notebook
    jupyter nbconvert --to notebook --execute --inplace confusionMatrxPlots.ipynb
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to execute confusionMatrxPlots.ipynb in iteration $iteration."
        exit 1
    fi

    # Calculate error metrics
    base_path="grassHmmDatabase/grassfam1"

    A=$(tail -n +2 "$base_path"/filtered_confusion_df | awk '$2 != $3' | wc -l)
    B=$(tail -n +2 "$base_path"/geneInfo | wc -l)
    C=$(tail -n +2 "$base_path"/filtered_confusion_df | wc -l)

    echo -e "mishits: $A"
    nohits=$((B - C))
    echo -e "nohits: $nohits"
    echo -e "total: $B"
    percent_error=$(awk -v A="$A" -v B="$B" -v C="$C" 'BEGIN {printf "%.8f", ((A + B - C) / B) * 100}')
    echo -e "%error=$percent_error%"

    if [[ -n "$prev_error" ]]; then
        diff=$(awk -v e1="$percent_error" -v e2="$prev_error" 'BEGIN {print (e1 > e2) ? e1 - e2 : e2 - e1}')
        threshold=0.000001
        if (( $(echo "$diff < $threshold" | bc -l) )); then
            echo "Error percentage change below threshold ($threshold). Stopping."

            cd ..
            ((iteration++))

            rm -rf "HmmRun$((iteration))"   
            while read -r item; do
                rm -f "${dir}/grassfamHmmLib/${item}.hmm"
            done < <(tail -n +2 "$dir/$base_path"/disagreeMatrix3 | awk '{print $1}')
            mkdir -p HMMFinal
            cat $dir/grassfamHmmLib/* > HMMFinal/grassfam
            hmmpress HMMFinal/grassfam

            # in confusionMatrix of dir, remove rows that contain items in disagreeMatrix3, then rerun confusionMatrxPlots.ipynb
            cd "$dir"
            awk 'NR==FNR {exclude[$1]; next} !($2 in exclude)' "$base_path"/disagreeMatrix3 "$base_path"/confusionMatrix > "$base_path"/temp
            mv "$base_path"/temp "$base_path"/confusionMatrix
            jupyter nbconvert --to notebook --execute --inplace confusionMatrxPlots.ipynb
            if [[ $? -ne 0 ]]; then
                echo "Error: Failed to execute confusionMatrxPlots.ipynb in final pass."
                exit 1
            fi
            cd ..

            rm -rf "HmmRun$((iteration))"
            break
        fi
    fi

    prev_error="$percent_error"

    # Go back to the parent directory to prepare for next iteration
    cd ..
    # Increment iteration
    ((iteration++))

    if (( iteration > 15 )); then
        echo "Reached maximum iteration limit (15). Stopping."

        rm -rf "HmmRun$((iteration))"
        while read -r item; do
            rm -f "${dir}/grassfamHmmLib/${item}.hmm"
        done < <(tail -n +2 "$dir/$base_path"/disagreeMatrix3 | awk '{print $1}')
        mkdir -p HMMFinal
        cat $dir/grassfamHmmLib/* > HMMFinal/grassfam
        hmmpress HMMFinal/grassfam

        # in confusionMatrix of dir, remove rows that contain items in disagreeMatrix3, then rerun confusionMatrxPlots.ipynb
        cd "$dir"
        awk 'NR==FNR {exclude[$1]; next} !($2 in exclude)' "$base_path"/disagreeMatrix3 "$base_path"/confusionMatrix > "$base_path"/temp
        mv "$base_path"/temp "$base_path"/confusionMatrix
        jupyter nbconvert --to notebook --execute --inplace confusionMatrxPlots.ipynb
        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to execute confusionMatrxPlots.ipynb in final pass."
            exit 1
        fi
        cd ..

        rm -rf "HmmRun$((iteration))"
        break
    fi
done

echo "GrassFam iterative clustering finished after $iteration iterations."