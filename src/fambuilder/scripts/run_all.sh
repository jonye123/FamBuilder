#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# run_all.sh
# Orchestrates the FamBuilder pipeline steps:
#   1. MAFFT alignment (maffter.sh)
#   2. HMM building (hmmbuilder.sh)
#   3. HMM pressing
#   4. Confusion matrix generation (confusionMatrixMaker.sh)
#
# Assumes it is called from inside HmmRun{N}/bashScripts/
# Part of the FamBuilder pipeline.
#
# Usage:
#   bash run_all.sh \
#       --project_name <name> \
#       --jobs <n> \
#       --conda_env <env_name>
# ---------------------------------------------------------------------------

# --- Argument parsing ---
project_name=""
jobs=8
conda_env=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project_name) project_name="$2"; shift 2 ;;
        --jobs)         jobs="$2";         shift 2 ;;
        --conda_env)    conda_env="$2";    shift 2 ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash run_all.sh --project_name <name> --jobs <n> --conda_env <env>"
            exit 1
            ;;
    esac
done

# --- Validate arguments ---
if [[ -z "$project_name" ]]; then
    echo "Error: --project_name is required."
    exit 1
fi

# --- Activate conda environment if specified ---
if [[ -n "$conda_env" ]]; then
    echo "Activating conda environment '$conda_env'..."

    # Support micromamba, mamba, and conda — use whichever is available
    if command -v micromamba &> /dev/null; then
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate "$conda_env"
    elif command -v mamba &> /dev/null; then
        eval "$(mamba shell hook --shell bash)"
        mamba activate "$conda_env"
    elif command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "$conda_env"
    else
        echo "Error: No conda/mamba/micromamba found in PATH."
        exit 1
    fi
fi

# --- Define paths relative to HmmRun{N}/ ---
# run_all.sh lives in HmmRun{N}/bashScripts/ so parent is HmmRun{N}/
run_dir="$(cd .. && pwd)"

hog_dir="${run_dir}/hogDirectoryFiltered"
msa_dir="${run_dir}/MSA"
hmm_lib_dir="${run_dir}/${project_name}HmmLib"
db_dir="${run_dir}/${project_name}HmmDatabase/${project_name}fam1"
hmm_lib="${run_dir}/${project_name}HmmDatabase/${project_name}"
scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log_dir="${scripts_dir}/bashLogs"

mkdir -p "$log_dir"
mkdir -p "$db_dir"

# --- Validate inputs ---
if [[ ! -d "$hog_dir" ]]; then
    echo "Error: hogDirectoryFiltered not found at '$hog_dir'."
    echo "Make sure hogFilter.py has been run before run_all.sh."
    exit 1
fi

echo "--- FamBuilder run_all.sh ---"
echo "Project     : $project_name"
echo "Jobs        : $jobs"
echo "Run dir     : $run_dir"
echo "HOG dir     : $hog_dir"
echo "MSA dir     : $msa_dir"
echo "HMM lib dir : $hmm_lib_dir"
echo "DB dir      : $db_dir"

# --- Step 1: MAFFT alignment ---
echo "deez nuts"
'''
echo "[1/4] Running MAFFT alignments..."
bash "${scripts_dir}/maffter.sh" \
    --input_dir "$hog_dir" \
    --output_dir "$msa_dir" \
    --jobs "$jobs" \
    ${prev_msa_dir:+--prev_msa_dir "$prev_msa_dir"} \
    ${prev_checksum_file:+--prev_checksum_file "$prev_checksum_file"} \
    > "$log_dir/maffter.log" 2> "$log_dir/maffter.err"

if [[ $? -ne 0 ]]; then
    echo "Error: maffter.sh failed. Check '$log_dir/maffter.err'."
    exit 1
fi
'''
# Verify MSA output exists
if [[ -z "$(ls -A "$msa_dir" 2>/dev/null)" ]]; then
    echo "Error: No MSA files produced in '$msa_dir'."
    exit 1
fi
echo "MAFFT alignments complete."

# --- Step 2: Build HMMs ---
echo ""
echo "[2/4] Building HMMs..."
bash "${scripts_dir}/hmmbuilder.sh" \
    --input_dir "$msa_dir" \
    --output_dir "$hmm_lib_dir" \
    --jobs "$jobs" \
    > "$log_dir/hmmbuilder.log" 2> "$log_dir/hmmbuilder.err"

if [[ $? -ne 0 ]]; then
    echo "Error: hmmbuilder.sh failed. Check '$log_dir/hmmbuilder.err'."
    exit 1
fi

# Verify HMMs were built
if [[ -z "$(ls -A "$hmm_lib_dir"/*.hmm 2>/dev/null)" ]]; then
    echo "Error: No HMM files produced in '$hmm_lib_dir'."
    exit 1
fi
echo "HMM building complete."

# --- Step 3: Press HMM library ---
echo ""
echo "[3/4] Pressing HMM library..."
hmm_count=$(ls "$hmm_lib_dir"/*.hmm | wc -l)
echo "Concatenating $hmm_count HMM files..."

cat "$hmm_lib_dir"/*.hmm > "$hmm_lib"

if [[ ! -s "$hmm_lib" ]]; then
    echo "Error: HMM library '$hmm_lib' is empty after concatenation."
    exit 1
fi

hmmpress "$hmm_lib" > "$log_dir/hmmpress.log" 2> "$log_dir/hmmpress.err"
if [[ $? -ne 0 ]]; then
    echo "Error: hmmpress failed. Check '$log_dir/hmmpress.err'."
    exit 1
fi
echo "HMM library pressed."

# --- Step 4: Confusion matrix ---
echo ""
echo "[4/4] Generating confusion matrix..."
bash "${scripts_dir}/confusionMatrixMaker.sh" \
    --hog_dir "$hog_dir" \
    --hmm_lib "$hmm_lib" \
    --output_dir "$db_dir" \
    --jobs "$jobs" \
    > "$log_dir/confusionMatrix.log" 2> "$log_dir/confusionMatrix.err"

if [[ $? -ne 0 ]]; then
    echo "Error: confusionMatrixMaker.sh failed. Check '$log_dir/confusionMatrix.err'."
    exit 1
fi

echo "Confusion matrix complete."
echo ""
echo "--- run_all.sh finished successfully ---"