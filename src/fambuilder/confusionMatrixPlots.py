#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description="Analyze HMM confusion matrix and prepare next iteration. "
                    "Part of the FamBuilder pipeline."
    )
    parser.add_argument("--project_name", required=True,
                        help="Name of the project (e.g. GrassFam, FungiFam). "
                             "Used to build output directory paths.")
    parser.add_argument("--iteration", type=int, required=True,
                        help="Current iteration number (1-based)")
    parser.add_argument("--iqr_multiplier", type=float, default=3.5,
                        help="IQR multiplier for score outlier filtering (default: 3.5)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Step 1: IQR filtering
# ---------------------------------------------------------------------------

def calculate_quartiles_and_filter(group, multiplier=3.5):
    """
    Filter out low-scoring genes within each HMM group using IQR.
    Genes scoring below Q1 - multiplier * IQR are removed.
    """
    Q1 = group["Score"].quantile(0.25)
    Q3 = group["Score"].quantile(0.75)
    IQR = Q3 - Q1
    median = group["Score"].median()
    lower_bound = Q1 - multiplier * IQR

    group = group.copy()
    group["Q1"] = Q1
    group["Q3"] = Q3
    group["IQR"] = IQR
    group["Median"] = median
    group["LowerBound"] = lower_bound

    return group[group["Score"] >= lower_bound]


def filter_confusion_matrix(confusion_matrix_path, output_path, iqr_multiplier):
    """
    Read raw confusionMatrix, apply IQR filtering per HMM group,
    and write filtered_confusion_df.
    Returns the filtered DataFrame.
    """
    print(f"Reading confusion matrix from '{confusion_matrix_path}'...")
    confusion_df = pd.read_csv(confusion_matrix_path, sep="\t")

    print(f"Applying IQR filtering (multiplier={iqr_multiplier})...")
    filtered_df = (
        confusion_df
        .groupby("HmmHogHit", group_keys=False)
        .apply(lambda g: calculate_quartiles_and_filter(g, iqr_multiplier))
        .reset_index(drop=True)
    )

    filtered_df.to_csv(output_path, sep="\t", index=False)
    print(f"Filtered confusion matrix written to '{output_path}' "
          f"({len(filtered_df)} / {len(confusion_df)} genes retained).")
    return filtered_df


# ---------------------------------------------------------------------------
# Step 2: Build disagree matrices
# ---------------------------------------------------------------------------

def build_disagree_matrices(filtered_confusion_df, db_path):
    """
    Build three disagreeMatrix files from the filtered confusion matrix.

    disagreeMatrix  — count of genes per (Hog, HmmHogHit) pair where they disagree
    disagreeMatrix2 — adds hogSize and ratio (mismatch / hogSize)
    disagreeMatrix3 — aggregates per HOG: total mismatch count, how many HMMs it
                      disagreed with, and overall ratio
    """
    print("Building disagree matrices...")

    # disagreeMatrix: pairs where Hog != HmmHogHit
    disagree_df = (
        filtered_confusion_df
        .groupby(["Hog", "HmmHogHit"])
        .size()
        .reset_index(name="count")
    )
    disagree_df = disagree_df[disagree_df["Hog"] != disagree_df["HmmHogHit"]]
    disagree_df = disagree_df[["count", "Hog", "HmmHogHit"]]
    disagree_df = disagree_df.rename(columns={
        "count": "count",
        "Hog": "hog",
        "HmmHogHit": "hmmHogHit"
    })
    disagree_df.to_csv(os.path.join(db_path, "disagreeMatrix"), sep="\t", index=False)

    # disagreeMatrix2: add hogSize and ratio
    hog_count = filtered_confusion_df["Hog"].value_counts()
    disagree_df["hogSize"] = disagree_df["hog"].map(hog_count)
    disagree_df["ratio"] = disagree_df["count"].values / disagree_df["hogSize"]
    disagree_df.rename(columns={"count": "mismatchCount"}, inplace=True)
    disagree_df = disagree_df[["hog", "hmmHogHit", "mismatchCount", "hogSize", "ratio"]]
    disagree_df.to_csv(os.path.join(db_path, "disagreeMatrix2"), sep="\t", index=False)

    # disagreeMatrix3: aggregate per HOG
    grouped = (
        disagree_df
        .groupby(["hog", "hogSize"])
        .agg({"mismatchCount": "sum"})
        .reset_index()
    )
    grouped["ratio"] = grouped["mismatchCount"] / grouped["hogSize"]
    grouped["hmmWrongCount"] = (
        disagree_df.groupby(["hog", "hogSize"]).size().values
    )
    grouped = grouped[["hog", "hmmWrongCount", "mismatchCount", "hogSize", "ratio"]]
    grouped = grouped.sort_values(by="mismatchCount", ascending=False)
    grouped.to_csv(os.path.join(db_path, "disagreeMatrix3"), sep="\t", index=False)

    print(f"Disagree matrices written to '{db_path}'.")
    return disagree_df, grouped


# ---------------------------------------------------------------------------
# Step 3: Build geneInfo
# ---------------------------------------------------------------------------

def build_gene_info(hog_directory, output_path):
    """
    Scan hogDirectoryFiltered and write a table of hog, gene, gene_length
    for every gene in every HOG FASTA file.
    """
    print(f"Building geneInfo from '{hog_directory}'...")

    with open(output_path, "w") as f_out:
        f_out.write("hog\tgene\tgene_length\n")

        for filename in os.listdir(hog_directory):
            if not filename.endswith(".fasta"):
                continue

            # Handle both N0.HOG001.fasta and HOG001.fasta formats
            hog_id = filename.replace(".fasta", "").replace("N0.", "")
            hog_file = os.path.join(hog_directory, filename)

            try:
                for record in SeqIO.parse(hog_file, "fasta"):
                    f_out.write(f"{hog_id}\t{record.id}\t{len(record.seq)}\n")
            except Exception as e:
                print(f"Warning: skipping corrupt file '{filename}': {e}")

    print(f"geneInfo written to '{output_path}'.")
    return pd.read_csv(output_path, sep="\t")


# ---------------------------------------------------------------------------
# Step 4: Build superMatrix
# ---------------------------------------------------------------------------

def build_super_matrix(filtered_confusion_df, gene_info_df, disagree_matrix3_df,
                        hog_msa_stats_path, output_path):
    """
    Join confusionMatrix, geneInfo, disagreeMatrix3, and hogMsaStats into
    a single superMatrix. Adds z-scores for gene length within each HOG.
    Also builds and returns wrongMatrix (rows where Hog != HmmHogHit).
    """
    print("Building superMatrix...")

    hog_msa_stats = pd.read_csv(hog_msa_stats_path, sep="\t")

    # Join all tables
    merged = pd.merge(filtered_confusion_df, gene_info_df,
                      left_on="Gene", right_on="gene")
    merged = pd.merge(merged, disagree_matrix3_df,
                      left_on="Hog", right_on="hog", how="left")
    merged = pd.merge(merged, hog_msa_stats,
                      left_on="Hog", right_on="HOG", how="left")

    # Drop redundant columns
    cols_to_drop = [c for c in ["gene", "hog_x", "hog_y", "HOG",
                                 "listOfLengths", "min", "max",
                                 "species_count", "gene_count"] if c in merged.columns]
    merged.drop(cols_to_drop, axis=1, inplace=True)

    # Calculate z-scores for gene length within each HOG
    gene_length = merged["gene_length"]
    mean = merged["mean"]
    sd = merged["sd"]

    # Replace infinite z-scores with 0 (happens when sd=0, i.e. all genes same length)
    z_scores = (gene_length - mean) / sd
    merged["z_scores"] = np.where(np.isfinite(z_scores), z_scores, 0)

    merged.to_csv(output_path, sep="\t", index=False)
    print(f"superMatrix written to '{output_path}' ({len(merged)} rows).")

    # wrongMatrix: genes hitting the wrong HMM
    wrong_matrix = merged[merged["Hog"] != merged["HmmHogHit"]]
    wrong_matrix_path = output_path.replace("superMatrix", "wrongMatrix")
    wrong_matrix.to_csv(wrong_matrix_path, sep="\t", index=False)
    print(f"wrongMatrix written to '{wrong_matrix_path}' ({len(wrong_matrix)} rows).")

    return merged, wrong_matrix


# ---------------------------------------------------------------------------
# Step 5: Find dropped genes (no HMM hit)
# ---------------------------------------------------------------------------

def find_dropped_genes(gene_info_df, filtered_confusion_df, output_path):
    """
    Find genes that are in geneInfo (i.e. they exist in hogDirectoryFiltered)
    but are absent from filtered_confusion_df (i.e. they scored to no HMM).
    These are 'no-hit' genes.
    """
    print("Finding no-hit (dropped) genes...")

    hog_set = set(gene_info_df["gene"])
    hmm_set = set(filtered_confusion_df["Gene"])
    no_hit_genes = hog_set - hmm_set

    dropped_df = gene_info_df[gene_info_df["gene"].isin(no_hit_genes)]
    dropped_df.to_csv(output_path, sep="\t", index=False)

    print(f"Found {len(no_hit_genes)} no-hit genes. Written to '{output_path}'.")
    return dropped_df



# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = get_args()

    # Build paths from project name
    db_path = f"{args.project_name}HmmDatabase/{args.project_name}fam1"
    hog_directory = "hogDirectoryFiltered"
    hog_msa_stats_path = "hogMsaStats.tsv"

    print(f"--- FamBuilder confusionMatrixPlots ---")
    print(f"Project  : {args.project_name}")
    print(f"Iteration: {args.iteration}")
    print(f"DB path  : {db_path}")

    # Validate inputs exist
    confusion_matrix_path = os.path.join(db_path, "confusionMatrix")
    if not os.path.exists(confusion_matrix_path):
        print(f"Error: confusionMatrix not found at '{confusion_matrix_path}'.")
        sys.exit(1)

    if not os.path.exists(hog_directory):
        print(f"Error: '{hog_directory}' not found.")
        sys.exit(1)

    if not os.path.exists(hog_msa_stats_path):
        print(f"Error: '{hog_msa_stats_path}' not found.")
        sys.exit(1)

    # Step 1: IQR filtering
    filtered_confusion_df = filter_confusion_matrix(
        confusion_matrix_path,
        os.path.join(db_path, "filtered_confusion_df"),
        args.iqr_multiplier
    )

    # Step 2: Disagree matrices
    disagree_df, disagree_matrix3_df = build_disagree_matrices(
        filtered_confusion_df, db_path
    )

    # Step 3: geneInfo
    gene_info_df = build_gene_info(
        hog_directory,
        os.path.join(db_path, "geneInfo")
    )

    # Step 4: superMatrix + wrongMatrix
    super_matrix_df, wrong_matrix_df = build_super_matrix(
        filtered_confusion_df,
        gene_info_df,
        disagree_matrix3_df,
        hog_msa_stats_path,
        os.path.join(db_path, "superMatrix")
    )

    # Step 5: Dropped genes
    find_dropped_genes(
        gene_info_df,
        filtered_confusion_df,
        os.path.join(db_path, "dropped")
    )


    print("confusionMatrixPlots complete.")


if __name__ == "__main__":
    main()