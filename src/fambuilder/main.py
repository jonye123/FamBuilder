#!/usr/bin/env python

import os
import sys
import json
import shutil
import argparse
import subprocess
import importlib.resources

import yaml


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description="FamBuilder: Iterative HMM-based gene family classification pipeline."
    )
    parser.add_argument("--config", required=True,
                        help="Path to config.yaml")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint if pipeline was interrupted")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(config_path):
    """Load and validate config.yaml."""
    if not os.path.exists(config_path):
        print(f"Error: Config file '{config_path}' not found.")
        sys.exit(1)

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Required fields
    required = ["project_name", "input", "fasta_dir"]
    for field in required:
        if field not in config:
            print(f"Error: Required field '{field}' missing from config.")
            sys.exit(1)

    # Validate input file exists
    if not os.path.exists(config["input"]):
        print(f"Error: Input file '{config['input']}' not found.")
        sys.exit(1)

    # Validate fasta_dir exists
    if not os.path.exists(config["fasta_dir"]):
        print(f"Error: FASTA directory '{config['fasta_dir']}' not found.")
        sys.exit(1)

    # Set defaults
    config.setdefault("min_species", 4)
    config.setdefault("min_groups", 0)
    config.setdefault("jobs", 8)
    config.setdefault("taxonomic_groups", {})
    config.setdefault("conda_env", "")
    config.setdefault("iqr_multiplier", 3.5)
    config.setdefault("max_iterations", 15)
    config.setdefault("convergence_threshold", 0.000001)
    config.setdefault("cleanup", True)

    return config


# ---------------------------------------------------------------------------
# Checkpointing
# ---------------------------------------------------------------------------

CHECKPOINT_FILE = "fambuilder_checkpoint.json"

def load_checkpoint():
    """Load checkpoint from disk if it exists."""
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r") as f:
            return json.load(f)
    return {"iteration": 0, "completed_steps": [], "prev_error": None}


def save_checkpoint(checkpoint):
    """Save checkpoint to disk."""
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(checkpoint, f, indent=2)


def mark_step_complete(checkpoint, step):
    """Mark a pipeline step as complete and save."""
    if step not in checkpoint["completed_steps"]:
        checkpoint["completed_steps"].append(step)
    save_checkpoint(checkpoint)


def step_is_complete(checkpoint, step):
    """Check if a step has already been completed."""
    return step in checkpoint["completed_steps"]


def reset_iteration_steps(checkpoint, iteration):
    """Clear completed steps for a given iteration (for fresh iteration start)."""
    checkpoint["completed_steps"] = [
        s for s in checkpoint["completed_steps"]
        if not s.startswith(f"iter{iteration}_")
    ]
    save_checkpoint(checkpoint)


# ---------------------------------------------------------------------------
# Script location
# ---------------------------------------------------------------------------

def get_scripts_dir():
    """
    Locate the bundled shell scripts directory from the installed package.
    Falls back to a local 'scripts/' directory for development use.
    """
    try:
        scripts_dir = importlib.resources.files("fambuilder") / "scripts"
        # Convert to a real path
        scripts_dir = str(scripts_dir)
        if os.path.isdir(scripts_dir):
            return scripts_dir
    except Exception:
        pass

    # Development fallback — look for scripts/ relative to this file
    local_scripts = os.path.join(os.path.dirname(__file__), "scripts")
    if os.path.isdir(local_scripts):
        return local_scripts

    print("Error: Could not locate bundled shell scripts.")
    sys.exit(1)


def get_python_scripts_dir():
    """
    Locate the bundled Python scripts from the installed package.
    Falls back to the directory containing this file for development use.
    """
    try:
        pkg_dir = importlib.resources.files("fambuilder")
        pkg_dir = str(pkg_dir)
        if os.path.isdir(pkg_dir):
            return pkg_dir
    except Exception:
        pass

    return os.path.dirname(__file__)


# ---------------------------------------------------------------------------
# Iteration setup
# ---------------------------------------------------------------------------

def setup_iteration_dir(iteration, config, scripts_dir, python_scripts_dir):
    """
    Create HmmRun{N}/ directory structure and copy all scripts into it.
    """
    run_dir = f"HmmRun{iteration}"
    bash_scripts_dir = os.path.join(run_dir, "bashScripts")
    bash_logs_dir = os.path.join(bash_scripts_dir, "bashLogs")

    os.makedirs(bash_scripts_dir, exist_ok=True)
    os.makedirs(bash_logs_dir, exist_ok=True)

    # Copy shell scripts
    for fname in os.listdir(scripts_dir):
        if fname.endswith(".sh"):
            shutil.copy(os.path.join(scripts_dir, fname), bash_scripts_dir)

    # Copy Python pipeline scripts
    for fname in ["hogFilter.py", "confusionMatrixPlots.py", "orthofinder_to_fasta.py"]:
        src = os.path.join(python_scripts_dir, fname)
        if os.path.exists(src):
            shutil.copy(src, run_dir)
        else:
            print(f"Warning: Could not find '{fname}' in '{python_scripts_dir}'.")

    # Copy input file
    shutil.copy(config["input"], run_dir)

    print(f"Iteration {iteration} directory ready at '{run_dir}'.")
    return run_dir


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def run_hog_filter(run_dir, config, iteration):
    """Run hogFilter.py for the current iteration."""
    cmd = [
        sys.executable,
        "hogFilter.py",
        "--input", os.path.basename(config["input"]),
        "--project_name", config["project_name"],
        "--iteration", str(iteration),
        "--min_species", str(config["min_species"]),
        "--min_groups", str(config["min_groups"]),
        "--taxonomic_groups", json.dumps(config["taxonomic_groups"]),
        "--jobs", str(config["jobs"]),
    ]

    if iteration == 1:
        cmd += ["--fasta_dir", config["fasta_dir"]]

    print(f"Running hogFilter.py (iteration {iteration})...")
    result = subprocess.run(cmd, cwd=run_dir)
    if result.returncode != 0:
        print(f"Error: hogFilter.py failed in iteration {iteration}.")
        sys.exit(1)


def run_all_sh(run_dir, config, iteration):
    """Run run_all.sh for the current iteration."""
    bash_scripts_dir = os.path.join(run_dir, "bashScripts")

    cmd = [
        "bash", "run_all.sh",
        "--project_name", config["project_name"],
        "--jobs", str(config["jobs"]),
    ]

    if config["conda_env"]:
        cmd += ["--conda_env", config["conda_env"]]

    # Pass previous MSA dir and checksum file for caching (iteration 2+)
    if iteration > 1:
        prev_run_dir = f"../HmmRun{iteration - 1}"
        prev_msa_dir = os.path.join(prev_run_dir, "MSA")
        prev_checksum_file = os.path.join(prev_msa_dir, "msa_checksums.tsv")

        # Use absolute paths to be safe
        prev_msa_abs = os.path.abspath(os.path.join(run_dir, prev_msa_dir))
        prev_checksum_abs = os.path.abspath(os.path.join(run_dir, prev_checksum_file))

        if os.path.isdir(prev_msa_abs):
            cmd += ["--prev_msa_dir", prev_msa_abs]
        if os.path.isfile(prev_checksum_abs):
            cmd += ["--prev_checksum_file", prev_checksum_abs]

    print(f"Running run_all.sh (iteration {iteration})...")
    result = subprocess.run(cmd, cwd=bash_scripts_dir)
    if result.returncode != 0:
        print(f"Error: run_all.sh failed in iteration {iteration}.")
        sys.exit(1)


def run_confusion_matrix_plots(run_dir, config, iteration):
    """Run confusionMatrixPlots.py for the current iteration."""
    cmd = [
        sys.executable,
        "confusionMatrixPlots.py",
        "--project_name", config["project_name"],
        "--iteration", str(iteration),
        "--iqr_multiplier", str(config["iqr_multiplier"]),
    ]

    print(f"Running confusionMatrixPlots.py (iteration {iteration})...")
    result = subprocess.run(cmd, cwd=run_dir)
    if result.returncode != 0:
        print(f"Error: confusionMatrixPlots.py failed in iteration {iteration}.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Error calculation
# ---------------------------------------------------------------------------

def calculate_error(run_dir, config):
    """
    Calculate % error for the current iteration.
    error = (mishits + nohits) / total_genes * 100

    mishits = genes where Hog != HmmHogHit in filtered_confusion_df
    nohits  = genes in geneInfo but not in filtered_confusion_df
    total   = all genes in geneInfo
    """
    db_path = os.path.join(
        run_dir,
        f"{config['project_name']}HmmDatabase",
        f"{config['project_name']}Outputs"
    )

    filtered_confusion_path = os.path.join(db_path, "filtered_confusion_df")
    gene_info_path = os.path.join(db_path, "geneInfo")

    if not os.path.exists(filtered_confusion_path) or not os.path.exists(gene_info_path):
        print("Error: Could not find filtered_confusion_df or geneInfo for error calculation.")
        sys.exit(1)

    # Count mishits (A): rows where Hog != HmmHogHit
    mishits = 0
    with open(filtered_confusion_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[1] != parts[2]:
                mishits += 1

    # Count total genes (B) and genes with HMM hits (C)
    total_genes = 0
    with open(gene_info_path) as f:
        next(f)  # skip header
        for _ in f:
            total_genes += 1

    hmm_hits = 0
    with open(filtered_confusion_path) as f:
        next(f)  # skip header
        for _ in f:
            hmm_hits += 1

    nohits = total_genes - hmm_hits

    if total_genes == 0:
        print("Error: No genes found in geneInfo.")
        sys.exit(1)

    percent_error = ((mishits + nohits) / total_genes) * 100

    print(f"  Mishits : {mishits}")
    print(f"  No-hits : {nohits}")
    print(f"  Total   : {total_genes}")
    print(f"  % Error : {percent_error:.8f}%")

    return percent_error


# ---------------------------------------------------------------------------
# Finalization
# ---------------------------------------------------------------------------

def finalize(run_dir, config, iteration):
    """
    Build the final HMM library:
    1. Remove the most disagreeing HOGs from the library
    2. Repress the library
    3. Run a final pass of confusionMatrixPlots.py
    4. Copy to HMMFinal/
    """
    print(f"\nFinalizing from '{run_dir}'...")

    project_name = config["project_name"]
    db_path = os.path.join(run_dir, f"{project_name}HmmDatabase", f"{project_name}Outputs")
    hmm_lib_dir = os.path.join(run_dir, f"{project_name}HmmLib")
    hmm_lib = os.path.join(run_dir, f"{project_name}HmmDatabase", project_name)
    disagree_matrix3 = os.path.join(db_path, "disagreeMatrix3")
    confusion_matrix = os.path.join(db_path, "confusionMatrix")

    # Remove most disagreeing HOMs from the HMM library
    if os.path.exists(disagree_matrix3):
        print("Removing most disagreeing HOGs from HMM library...")
        with open(disagree_matrix3) as f:
            next(f)  # skip header
            for line in f:
                hog_id = line.strip().split("\t")[0]
                hmm_file = os.path.join(hmm_lib_dir, f"{hog_id}.hmm")
                if os.path.exists(hmm_file):
                    os.remove(hmm_file)

    # Repress the library
    print("Rebuilding and pressing final HMM library...")
    hmm_files = [f for f in os.listdir(hmm_lib_dir) if f.endswith(".hmm")]
    if not hmm_files:
        print("Error: No HMM files remaining after removing disagreeing HOGs.")
        sys.exit(1)

    with open(hmm_lib, "w") as out_f:
        for hmm_file in hmm_files:
            with open(os.path.join(hmm_lib_dir, hmm_file)) as in_f:
                out_f.write(in_f.read())

    # Remove existing hmmpress index files if they exist
    for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
        index_file = hmm_lib + ext
        if os.path.exists(index_file):
            os.remove(index_file)

    result = subprocess.run(["hmmpress", hmm_lib])
    if result.returncode != 0:
        print("Error: hmmpress failed during finalization.")
        sys.exit(1)

    # Remove disagreeing HOGs from confusion matrix and rerun analysis
    print("Running final confusionMatrixPlots pass...")
    if os.path.exists(disagree_matrix3) and os.path.exists(confusion_matrix):
        # Load HOGs to exclude
        exclude_hogs = set()
        with open(disagree_matrix3) as f:
            next(f)
            for line in f:
                exclude_hogs.add(line.strip().split("\t")[0])

        # Filter confusion matrix
        temp_path = confusion_matrix + ".tmp"
        with open(confusion_matrix) as in_f, open(temp_path, "w") as out_f:
            header = next(in_f)
            out_f.write(header)
            for line in in_f:
                hog = line.strip().split("\t")[1]
                if hog not in exclude_hogs:
                    out_f.write(line)
        os.replace(temp_path, confusion_matrix)

    run_confusion_matrix_plots(run_dir, config, iteration)

    # Copy to HMMFinal/
    os.makedirs("HMMFinal", exist_ok=True)
    for ext in ["", ".h3f", ".h3i", ".h3m", ".h3p"]:
        src = hmm_lib + ext
        if os.path.exists(src):
            shutil.copy(src, "HMMFinal/")

    print(f"Final HMM library written to 'HMMFinal/'.")


def cleanup_intermediate_dirs(final_iteration):
    """Delete all HmmRun{N} directories after successful completion."""
    print("Cleaning up intermediate directories...")
    for i in range(1, final_iteration):
        run_dir = f"HmmRun{i}"
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)
            print(f"  Removed '{run_dir}'.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = get_args()
    config = load_config(args.config)

    print("--- FamBuilder ---")
    print(f"Project   : {config['project_name']}")
    print(f"Input     : {config['input']}")
    print(f"FASTA dir : {config['fasta_dir']}")
    print(f"Jobs      : {config['jobs']}")
    print(f"Max iter  : {config['max_iterations']}")
    print(f"Threshold : {config['convergence_threshold']}")

    # Locate bundled scripts
    scripts_dir = get_scripts_dir()
    python_scripts_dir = get_python_scripts_dir()

    # Load or initialize checkpoint
    if args.resume:
        checkpoint = load_checkpoint()
        iteration = checkpoint["iteration"] + 1
        reset_iteration_steps(checkpoint, iteration)
    else:
        checkpoint = {"iteration": 0, "completed_steps": [], "prev_error": None}
        iteration = 1
    save_checkpoint(checkpoint)
    prev_error = checkpoint["prev_error"]

    while True:
        print(f"\n{'='*50}")
        print(f"Iteration {iteration}")
        print(f"{'='*50}")

        run_dir = f"HmmRun{iteration}"

        # Step 1: Setup iteration directory
        step_setup = f"iter{iteration}_setup"
        if not step_is_complete(checkpoint, step_setup):
            setup_iteration_dir(iteration, config, scripts_dir, python_scripts_dir)
            mark_step_complete(checkpoint, step_setup)
        else:
            print(f"Skipping setup (already complete).")

        # Step 2: hogFilter
        step_hog = f"iter{iteration}_hogFilter"
        if not step_is_complete(checkpoint, step_hog):
            run_hog_filter(run_dir, config, iteration)
            mark_step_complete(checkpoint, step_hog)
        else:
            print(f"Skipping hogFilter (already complete).")

        # Step 3: run_all.sh
        step_run_all = f"iter{iteration}_run_all"
        if not step_is_complete(checkpoint, step_run_all):
            run_all_sh(run_dir, config, iteration)
            mark_step_complete(checkpoint, step_run_all)
        else:
            print(f"Skipping run_all.sh (already complete).")

        # Step 4: confusionMatrixPlots
        step_confusion = f"iter{iteration}_confusionMatrixPlots"
        if not step_is_complete(checkpoint, step_confusion):
            run_confusion_matrix_plots(run_dir, config, iteration)
            mark_step_complete(checkpoint, step_confusion)
        else:
            print(f"Skipping confusionMatrixPlots (already complete).")

        # Step 5: Calculate error
        print(f"\nCalculating error for iteration {iteration}...")
        percent_error = calculate_error(run_dir, config)

        # Step 6: Check convergence
        converged = False
        if prev_error is not None:
            diff = abs(percent_error - prev_error)
            print(f"  Error change: {diff:.8f} (threshold: {config['convergence_threshold']})")
            if diff < config["convergence_threshold"]:
                print(f"Converged after {iteration} iterations.")
                converged = True

        if iteration >= config["max_iterations"]:
            print(f"Reached maximum iteration limit ({config['max_iterations']}).")
            converged = True

        if converged:
            # Finalize
            finalize(run_dir, config, iteration)

            # Cleanup intermediate dirs if requested
            if config["cleanup"]:
                cleanup_intermediate_dirs(iteration)

            # Remove checkpoint on successful completion
            if os.path.exists(CHECKPOINT_FILE):
                os.remove(CHECKPOINT_FILE)

            print(f"\nFamBuilder finished successfully after {iteration} iterations.")
            print(f"Final HMM library is in 'HMMFinal/'.")
            break

        # Prepare for next iteration
        prev_error = percent_error
        checkpoint["prev_error"] = percent_error
        checkpoint["iteration"] = iteration
        save_checkpoint(checkpoint)

        iteration += 1


if __name__ == "__main__":
    main()