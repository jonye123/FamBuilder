#!/usr/bin/env python

import os
import sys
import argparse
import shutil
import statistics

import pandas as pd
from Bio import SeqIO

try:
    import orthofinder_to_fasta
except ImportError:
    print("Warning: orthofinder_to_fasta.py not found. Iteration 1 extraction will fail.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description="Filter HOGs and generate stats. Part of the FamBuilder pipeline.")
    parser.add_argument("--input", required=True,
                        help="Path to N0.tsv (OrthoFinder HOG file)")
    parser.add_argument("--fasta_dir",
                        help="Directory containing per-species proteome FASTA files "
                             "(required for iteration 1)")
    parser.add_argument("--project_name", required=True,
                        help="Name of the project (e.g. GrassFam, FungiFam). "
                             "Used to name output directories.")
    parser.add_argument("--iteration", type=int, required=True,
                        help="Current iteration number (1-based)")
    parser.add_argument("--min_species", type=int, default=4,
                        help="Minimum number of species required to keep a HOG (default: 4)")
    parser.add_argument("--min_groups", type=int, default=0,
                        help="Minimum number of taxonomic groups represented in a HOG. "
                             "Used together with --min_species for group-aware filtering. "
                             "(default: 0, disabled)")
    parser.add_argument("--taxonomic_groups", type=str, default="{}",
                        help="JSON string defining taxonomic groups." 
                            "Only used when --min_groups > 0.")
    parser.add_argument("--jobs", type=int, default=8,
                        help="(Informational) Number of parallel jobs — passed through "
                             "to shell scripts by main.sh")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

def calculate_stats(hog_directory, output_tsv, species_list):
    """
    Scan a directory of FASTA files and write a statistics TSV.
    Returns the stats as a DataFrame, or an empty DataFrame if no files found.

    species_list is passed in explicitly so that the column set is consistent
    across calls (avoids the 'sample 100 files' bug from the earlier draft).
    """
    if not os.path.exists(hog_directory):
        print(f"Directory '{hog_directory}' does not exist — writing empty stats.")
        return pd.DataFrame()

    fasta_files = [f for f in os.listdir(hog_directory) if f.endswith(".fasta")]
    if not fasta_files:
        print(f"No FASTA files found in '{hog_directory}'.")
        return pd.DataFrame()

    rows = []

    for filename in fasta_files:
        hog_id = filename.replace(".fasta", "")
        hog_file = os.path.join(hog_directory, filename)

        species_found = set()
        lengths = []

        try:
            for record in SeqIO.parse(hog_file, "fasta"):
                lengths.append(len(record.seq))
                parts = record.id.split("|")
                if parts:
                    species_found.add(parts[0])
        except Exception as e:
            print(f"Skipping corrupt file '{filename}': {e}")
            continue

        if not lengths:
            continue

        gene_count = len(lengths)
        sd_len = statistics.stdev(lengths) if gene_count > 1 else 0

        species_binary = [1 if s in species_found else 0 for s in species_list]

        rows.append({
            "HOG": hog_id,
            "min": min(lengths),
            "max": max(lengths),
            "median": statistics.median(lengths),
            "mean": statistics.mean(lengths),
            "sd": sd_len,
            "species_count": len(species_found),
            "gene_count": gene_count,
            "listOfLengths": ",".join(map(str, lengths)),
            **{s: b for s, b in zip(species_list, species_binary)},
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Stats written to '{output_tsv}' ({len(df)} HOGs).")
    return df


# ---------------------------------------------------------------------------
# Species list detection
# ---------------------------------------------------------------------------

def detect_species_list(hog_directory):
    """
    Scan ALL FASTA files in hog_directory to build a complete, sorted species list.
    This is done once upfront so every downstream function uses the same set.
    """
    all_species = set()
    for filename in os.listdir(hog_directory):
        if not filename.endswith(".fasta"):
            continue
        try:
            for record in SeqIO.parse(os.path.join(hog_directory, filename), "fasta"):
                parts = record.id.split("|")
                if parts:
                    all_species.add(parts[0])
        except Exception:
            pass
    return sorted(all_species)


# ---------------------------------------------------------------------------
# HOG splitting (restored from original notebook logic)
# ---------------------------------------------------------------------------

def split_and_clean_hogs(hog_directory, gene_and_hog_to_remove, zeros_to_remove):
    """
    For each FASTA file in hog_directory:

    1. Remove genes that were 'no-hits' (present in geneInfo but absent from
       filtered_confusion_df). These are stored in zeros_to_remove keyed by hog_id.

    2. For genes that hit the WRONG HMM, split them out into a new sub-family file
       rather than discarding them. gene_and_hog_to_remove is a dict keyed by
       (source_hog_id, wrong_hmm_id) -> [list of gene IDs].

       Example: if HOG001 has 3 genes that hit HOG002's HMM, those 3 genes are
       written to HOG001-1.fasta, and the remaining genes stay in HOG001-0.fasta.

    3. Clean up empty -0 files: if the original HOG is entirely emptied (ratio=1.0),
       delete the empty -0 file and rename the -1 file back to the base name.
    """
    print("Splitting and cleaning HOG FASTA files...")

    for filename in os.listdir(hog_directory):
        if not filename.endswith(".fasta"):
            continue

        filepath = os.path.join(hog_directory, filename)

        # Parse the hog_id — handle both N0.HOG001.fasta and HOG001.fasta formats
        # Also handle suffixed files like N0.HOG001-1.fasta
        basename = filename.replace(".fasta", "")
        # Extract the core HOG id (strip N0. prefix if present, then take part before -)
        hog_id_full = basename.replace("N0.", "")
        hog_id = hog_id_full.split("-")[0]

        # Load all sequences into memory
        sequences = {}
        try:
            for record in SeqIO.parse(filepath, "fasta"):
                sequences[record.id] = str(record.seq)
        except Exception as e:
            print(f"Skipping corrupt file '{filename}': {e}")
            continue

        # Step 1: Remove no-hit genes (zeros_to_remove)
        for gene in zeros_to_remove.get(hog_id, []):
            sequences.pop(gene, None)

        # Step 2: Split out disagreeing genes into sub-family files
        file_counter = 0
        for (hog_key1, hog_key2), genes in gene_and_hog_to_remove.items():
            if hog_id != hog_key1:
                continue

            file_counter += 1
            removed_sequences = {}
            for gene in genes:
                if gene in sequences:
                    removed_sequences[gene] = sequences.pop(gene)
                else:
                    print(f"Warning: gene '{gene}' listed for removal not found in '{filename}'")

            if removed_sequences:
                new_path = os.path.join(
                    hog_directory,
                    f"{os.path.splitext(filename)[0]}-{file_counter}.fasta"
                )
                with open(new_path, "w") as f:
                    for gene_id, seq in removed_sequences.items():
                        f.write(f">{gene_id}\n{seq}\n")

        # Step 3: Write back the (possibly reduced) original sequences
        if file_counter == 0:
            # No splits were made — overwrite in place with the same filename
            with open(filepath, "w") as f:
                for gene_id, seq in sequences.items():
                    f.write(f">{gene_id}\n{seq}\n")
        else:
            # Splits were made — write original remainder as -0
            original_new_path = os.path.join(
                hog_directory,
                f"{os.path.splitext(filename)[0]}-0.fasta"
            )
            with open(original_new_path, "w") as f:
                for gene_id, seq in sequences.items():
                    f.write(f">{gene_id}\n{seq}\n")
            # Remove the unsuffixed original
            os.remove(filepath)

    # Step 4: Clean up empty -0 files (these are HOGs where ratio=1.0,
    # meaning ALL genes were reassigned — rename the -1 back to the base name)
    for filename in os.listdir(hog_directory):
        if not filename.endswith("-0.fasta"):
            continue

        filepath = os.path.join(hog_directory, filename)
        hog_id = filename.split(".")[1].split("-")[0] if "." in filename else filename.split("-")[0]

        if os.path.getsize(filepath) == 0:
            print(f"Empty -0 file detected, cleaning up: {filename}")
            os.remove(filepath)

            # Rename the corresponding -1 file back to the base name
            for related_file in os.listdir(hog_directory):
                related_hog_id = (
                    related_file.split(".")[1].split("-")[0]
                    if "." in related_file
                    else related_file.split("-")[0]
                )
                if related_hog_id == hog_id and related_file.endswith("-1.fasta"):
                    new_name = related_file.replace("-1.fasta", ".fasta")
                    os.rename(
                        os.path.join(hog_directory, related_file),
                        os.path.join(hog_directory, new_name),
                    )
                    break

    print("HOG splitting and cleaning complete.")


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------
# Taxonomic groups for group-aware filtering.
# This is populated at runtime from --taxonomic_groups (set via config.yaml).
TAXONOMIC_GROUPS = {}

def count_groups_represented(row):
    """
    Count how many taxonomic groups have at least one species present in the HOG.
    If TAXONOMIC_GROUPS is empty, falls back to counting individual species present,
    so min_groups behaves like min_species in that case.
    """
    if TAXONOMIC_GROUPS:
        return sum(
            1 for members in TAXONOMIC_GROUPS.values()
            if any(int(row.get(s, 0)) for s in members)
        )
    else:
        return int(row.get("species_count", 0))

def passes_filter(row, min_species, min_groups=0):
    """Return True if a HOG row should be kept."""
    species_ok = row["species_count"] >= min_species
    
    if min_groups > 0:
        groups_ok = count_groups_represented(row) >= min_groups
        return species_ok and groups_ok
    else:
        return species_ok

def copy_filtered_hogs(hog_stats_df, hog_directory, output_dir, min_species, min_groups=0):
    """
    Copy HOG FASTA files that pass the filter into output_dir.
    Handles both N0-prefixed and plain filenames, and suffixed files (HOG001-1.fasta).
    """
    os.makedirs(output_dir, exist_ok=True)
    count_kept = 0

    for _, row in hog_stats_df.iterrows():
        if not passes_filter(row, min_species, min_groups):
            continue

        hog_name = str(row["HOG"])

        # Try the filename as stored in the stats (already has N0. stripped by calculate_stats)
        candidates = [
            os.path.join(hog_directory, f"{hog_name}.fasta"),
            os.path.join(hog_directory, f"N0.{hog_name}.fasta"),
        ]

        copied = False
        for src in candidates:
            if os.path.exists(src):
                shutil.copy(src, output_dir)
                count_kept += 1
                copied = True
                break

        if not copied:
            print(f"Warning: FASTA file not found for HOG '{hog_name}'")

    print(f"Filtering complete — kept {count_kept} / {len(hog_stats_df)} HOGs.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = get_args()

    # Validate fasta_dir is provided for iteration 1
    if args.iteration == 1 and not args.fasta_dir:
        print("Error: --fasta_dir is required for iteration 1.")
        sys.exit(1)

    # Build database path using project name
    db_path = f"{args.project_name}HmmDatabase/{args.project_name}Outputs"

    print(f"--- FamBuilder hogFilter ---")
    print(f"Project : {args.project_name}")
    print(f"Iteration: {args.iteration}")
    print(f"DB path  : {db_path}")

    # ------------------------------------------------------------------
    # ITERATION 1: Extract HOG FASTAs from N0.tsv
    # ------------------------------------------------------------------
    if args.iteration == 1:
        print("Iteration 1: extracting sequences from OrthoFinder output...")

        if not os.path.exists(args.input):
            print(f"Error: Input file '{args.input}' not found.")
            sys.exit(1)

        try:
            orthogroups = orthofinder_to_fasta.parse_orthogroups(args.input)
            orthofinder_to_fasta.extract_sequences(orthogroups, args.fasta_dir, "./hogDirectory")
            print("HOG FASTAs written to ./hogDirectory")
        except AttributeError:
            print("Error: orthofinder_to_fasta module is missing required functions.")
            sys.exit(1)

    # ------------------------------------------------------------------
    # ITERATION 2+: Refine HOGs from previous run
    # ------------------------------------------------------------------
    else:
        print(f"Iteration {args.iteration}: refining HOGs from previous run...")

        prev_run = f"../HmmRun{args.iteration - 1}"
        # Use the same db path pattern for the previous run
        prev_db = f"{prev_run}/{args.project_name}HmmDatabase/{args.project_name}fam1"

        if not os.path.exists(f"{prev_run}/hogDirectoryFiltered"):
            print(f"Error: Previous filtered directory not found at '{prev_run}/hogDirectoryFiltered'.")
            sys.exit(1)

        # Start fresh from the previous run's filtered directory
        if os.path.exists("./hogDirectory"):
            shutil.rmtree("./hogDirectory")
        shutil.copytree(f"{prev_run}/hogDirectoryFiltered", "./hogDirectory")
        print("Copied filtered directory from previous run.")

        # --- Build gene_and_hog_to_remove from superMatrix ---
        # These are genes where Hog != HmmHogHit (they hit the wrong HMM)
        gene_and_hog_to_remove = {}
        sm_path = os.path.join(prev_db, "superMatrix")
        if os.path.exists(sm_path):
            sm = pd.read_csv(sm_path, sep="\t")
            disagreements = sm[sm["Hog"] != sm["HmmHogHit"]]
            gene_and_hog_to_remove = (
                disagreements.groupby(["Hog", "HmmHogHit"])["Gene"]
                .apply(list)
                .to_dict()
            )
            print(f"Found {len(disagreements)} disagreement genes across "
                  f"{len(gene_and_hog_to_remove)} (Hog, HmmHit) pairs.")
        else:
            print(f"Warning: superMatrix not found at '{sm_path}'. Skipping disagreement removal.")

        # --- Build zeros_to_remove from geneInfo vs filtered_confusion_df ---
        # These are genes that scored to nothing (no HMM hit at all)
        zeros_to_remove = {}
        gene_info_path = os.path.join(prev_db, "geneInfo")
        confusion_path = os.path.join(prev_db, "filtered_confusion_df")
        if os.path.exists(gene_info_path) and os.path.exists(confusion_path):
            gene_info = pd.read_csv(gene_info_path, sep="\t")
            hmm_list = set(pd.read_csv(confusion_path, sep="\t")["Gene"])
            hog_list = set(gene_info["gene"])
            no_hit_genes = hog_list - hmm_list
            filtered_gene_info = gene_info[gene_info["gene"].isin(no_hit_genes)]
            zeros_to_remove = (
                filtered_gene_info.groupby("hog")["gene"].apply(list).to_dict()
            )
            print(f"Found {len(no_hit_genes)} no-hit genes to remove.")
        else:
            print("Warning: geneInfo or filtered_confusion_df not found. Skipping no-hit removal.")

        # --- Apply splitting and cleaning ---
        split_and_clean_hogs("./hogDirectory", gene_and_hog_to_remove, zeros_to_remove)

    # ------------------------------------------------------------------
    # Detect species list from the full hogDirectory
    # ------------------------------------------------------------------
    print("Detecting species list from hogDirectory...")
    species_list = detect_species_list("./hogDirectory")
    if not species_list:
        print("Error: No species detected from hogDirectory. Check FASTA sequence ID format (expected: SpeciesCode|GeneID).")
        sys.exit(1)
    print(f"Detected {len(species_list)} species: {species_list}")

    # ------------------------------------------------------------------
    # Calculate stats
    # ------------------------------------------------------------------
    print("Generating hogMsaStats.tsv...")
    hog_stats_df = calculate_stats("./hogDirectory", "hogMsaStats.tsv", species_list)

    if hog_stats_df.empty:
        print("Error: No HOGs found to process.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Filter HOGs
    # ------------------------------------------------------------------
    if args.min_groups > 0:
        filter_mode = f"group-aware (min_species >= {args.min_species} && min_groups >= {args.min_groups})"
    else:
        filter_mode = f"generic (min_species >= {args.min_species})"
    print(f"Filtering HOGs using {filter_mode} filter...")
    copy_filtered_hogs(hog_stats_df, "./hogDirectory", "./hogDirectoryFiltered",
                       args.min_species, args.min_groups)

    # ------------------------------------------------------------------
    # Recalculate stats for filtered set
    # ------------------------------------------------------------------
    print("Generating hogMsaStats2.tsv for filtered set...")
    calculate_stats("./hogDirectoryFiltered", "hogMsaStats2.tsv", species_list)

    print("hogFilter complete.")


if __name__ == "__main__":
    main()