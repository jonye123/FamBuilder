#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import norm
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description="FamBuilder: Visualize HOG filtering and HMM confusion results."
    )
    parser.add_argument("--project_name", required=True,
                        help="Name of the project (e.g. GrassFam, FungiFam)")
    parser.add_argument("--iteration", type=int, required=True,
                        help="Iteration number to visualize")
    parser.add_argument("--plot", default="all",
                        choices=["all", "hog_stats", "confusion", "network"],
                        help="Which plots to generate (default: all)")
    parser.add_argument("--output_dir", default=None,
                        help="Directory to save plots (default: show interactively)")
    parser.add_argument("--iqr_multiplier", type=float, default=3.5,
                        help="IQR multiplier used in filtering (default: 3.5)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def save_or_show(fig, output_dir, filename):
    """Save figure to output_dir or show interactively."""
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, filename)
        fig.savefig(path, bbox_inches="tight", dpi=150)
        print(f"Saved: {path}")
        plt.close(fig)
    else:
        plt.show()


def load_hog_stats(run_dir):
    path = os.path.join(run_dir, "hogMsaStats2.tsv")
    if not os.path.exists(path):
        print(f"Warning: '{path}' not found. Skipping HOG stats plots.")
        return None
    return pd.read_csv(path, sep="\t")


def load_confusion_data(run_dir, db_path):
    """Load all confusion matrix related files."""
    data = {}

    for name, filename in [
        ("confusion_df", "confusionMatrix"),
        ("filtered_confusion_df", "filtered_confusion_df"),
        ("disagree_df", "disagreeMatrix2"),
        ("super_matrix", "superMatrix"),
        ("wrong_matrix", "wrongMatrix"),
    ]:
        path = os.path.join(db_path, filename)
        if os.path.exists(path):
            data[name] = pd.read_csv(path, sep="\t")
        else:
            print(f"Warning: '{path}' not found.")
            data[name] = None

    return data


# ---------------------------------------------------------------------------
# HOG Stats Plots (from hogFilter.ipynb)
# ---------------------------------------------------------------------------

def plot_hog_stats(hog_stats, output_dir, grass_filter=False):
    """Generate all HOG stats plots."""
    print("Generating HOG stats plots...")

    pacmad = ["Svi", "Pvi", "Pha", "Msi", "Sbi", "Zma", "Ecu", "Pau"]
    bop = ["Osa", "Ogl", "Ola", "Lpe", "Hvu", "Ata", "Tae", "Bdi"]

    # 1. Line plot: number of HOGs with at least X species
    fig, ax = plt.subplots()
    list_xy = [(i, len(hog_stats[hog_stats["species_count"] >= i]))
               for i in range(0, hog_stats["species_count"].max() + 2)]
    ax.plot(*zip(*list_xy))
    ax.set_xlabel("Minimum species count (inclusive)")
    ax.set_ylabel("Number of HOGs")
    ax.set_title("Number of HOGs with at least X species")
    ax.set_ylim(0)
    ax.grid()
    save_or_show(fig, output_dir, "hog_species_cutoff.png")

    # 2. Histogram: species count distribution
    fig, ax = plt.subplots()
    x_min, x_max = 0, int(hog_stats["species_count"].max())
    filtered = hog_stats["species_count"][
        (hog_stats["species_count"] >= x_min) &
        (hog_stats["species_count"] <= x_max)
    ]
    ax.hist(filtered, bins=x_max - x_min, range=(x_min, x_max))
    ax.set_title("Histogram of species count in HOG")
    ax.set_xlabel("Species count")
    ax.set_ylabel("Frequency")
    ax.set_xlim(x_min, x_max)
    ax.grid()
    save_or_show(fig, output_dir, "hog_species_histogram.png")

    # 3. Histogram: gene count distribution
    fig, ax = plt.subplots()
    x_min, x_max = 0, min(50, int(hog_stats["gene_count"].max()))
    filtered = hog_stats["gene_count"][
        (hog_stats["gene_count"] >= x_min) &
        (hog_stats["gene_count"] <= x_max)
    ]
    ax.hist(filtered, bins=x_max - x_min, range=(x_min, x_max))
    ax.set_title("Histogram of gene count in HOG")
    ax.set_xlabel("Gene count")
    ax.set_ylabel("Frequency")
    ax.set_xlim(x_min, x_max)
    ax.grid()
    save_or_show(fig, output_dir, "hog_gene_histogram.png")

    # 4 & 5. BOP/PACMAD plots (only if grass species columns are present)
    if grass_filter and all(s in hog_stats.columns for s in pacmad + bop):
        hog_stats = hog_stats.copy()
        hog_stats["BOP"] = hog_stats[bop].sum(axis=1)
        hog_stats["PACMAD"] = hog_stats[pacmad].sum(axis=1)

        frequency_table = (
            hog_stats.groupby(["BOP", "PACMAD"])
            .size()
            .unstack(fill_value=0)
            .reindex(index=range(9), columns=range(9), fill_value=0)
            .sort_index(axis=0, ascending=False)
            .sort_index(axis=1)
        )

        # 4. 2D heatmap
        fig, ax = plt.subplots()
        sns.heatmap(frequency_table, annot=True, fmt="d",
                    cbar_kws={"label": "Frequency"}, cmap="viridis", ax=ax)
        ax.set_title("Heatmap of BOP vs PACMAD Species Represented in HOGs")
        ax.set_xlabel("PACMAD Count")
        ax.set_ylabel("BOP Count")
        plt.tight_layout()
        save_or_show(fig, output_dir, "bop_pacmad_heatmap.png")

        # 5. 3D bar chart
        fig = plt.figure(figsize=(10, 7))
        ax3d = fig.add_subplot(111, projection="3d")
        x, y = np.meshgrid(frequency_table.columns, frequency_table.index)
        x, y = x.ravel(), y.ravel()
        z = np.zeros_like(x)
        dz = frequency_table.values.ravel()
        ax3d.bar3d(x, y, z, dx=1, dy=1, dz=dz, shade=True)
        ax3d.set_xlabel("PACMAD Count")
        ax3d.set_ylabel("BOP Count")
        ax3d.set_zlabel("Frequency")
        ax3d.set_title("3D Histogram of BOP vs PACMAD Species in HOGs")
        plt.tight_layout()
        save_or_show(fig, output_dir, "bop_pacmad_3d.png")


# ---------------------------------------------------------------------------
# Confusion Matrix Plots (from confusionMatrxPlots.ipynb)
# ---------------------------------------------------------------------------

def plot_confusion(data, output_dir, iqr_multiplier):
    """Generate all confusion matrix related plots."""
    print("Generating confusion matrix plots...")

    confusion_df = data.get("confusion_df")
    filtered_df = data.get("filtered_confusion_df")
    disagree_df = data.get("disagree_df")
    super_matrix = data.get("super_matrix")
    wrong_matrix = data.get("wrong_matrix")

    # 6. Histogram: all HMM scores
    if confusion_df is not None:
        fig, ax = plt.subplots()
        scores = confusion_df["Score"]
        num_bins = max(1, int(scores.max()) - int(scores.min()) // 2 + 1)
        ax.hist(scores, bins=num_bins)
        ax.set_title("Histogram of HMM Scores")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency")
        ax.grid()
        save_or_show(fig, output_dir, "score_histogram_all.png")

        # 7. Histogram: erroneous hits only
        erroneous = confusion_df[confusion_df["Hog"] != confusion_df["HmmHogHit"]]
        fig, ax = plt.subplots()
        num_bins = max(1, int(erroneous["Score"].max()) - int(erroneous["Score"].min()) // 2 + 1)
        ax.hist(erroneous["Score"], bins=num_bins)
        ax.set_title("Score of Erroneous HMM Hits")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency")
        ax.grid()
        save_or_show(fig, output_dir, "score_histogram_erroneous.png")

        # 8. Histogram: consensus vs erroneous overlaid
        consensus = confusion_df[confusion_df["Hog"] == confusion_df["HmmHogHit"]]
        fig, ax = plt.subplots()
        num_bins = max(1, int(confusion_df["Score"].max()) - int(confusion_df["Score"].min()) // 2 + 1)
        ax.hist(consensus["Score"], bins=num_bins, label="Consensus Hit")
        ax.hist(erroneous["Score"], bins=num_bins, color="red", label="Erroneous Hit", alpha=0.7)
        ax.set_title("HMM Score: Consensus vs Erroneous Hits")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency")
        ax.legend()
        ax.grid()
        save_or_show(fig, output_dir, "score_histogram_overlaid.png")

    # 9. Scatter: IQR vs median score per HMM
    if filtered_df is not None:
        median_score = np.log10(filtered_df.groupby("HmmHogHit")["Score"].median())
        bottom_q = filtered_df.groupby("HmmHogHit")["Score"].quantile(0.25)
        top_q = filtered_df.groupby("HmmHogHit")["Score"].quantile(0.75)
        iqr = np.log10(top_q - bottom_q + 1)

        m, b = np.polyfit(median_score, iqr, deg=1)

        fig, ax = plt.subplots()
        ax.scatter(median_score, iqr, s=0.1, alpha=0.3)
        ax.plot(median_score, m * median_score + b, label=f"y={m:.2f}x+{b:.2f}")
        ax.set_title("IQR vs Median Score per HMM")
        ax.set_xlabel("Log Median Score")
        ax.set_ylabel("Log IQR")
        ax.legend()
        ax.grid()
        save_or_show(fig, output_dir, "iqr_vs_median_score.png")

    # 10 & 11. Scatter plots using superMatrix and wrongMatrix
    if super_matrix is not None and wrong_matrix is not None:

        # 10. Median vs mean gene length
        fig, ax = plt.subplots()
        ax.scatter(super_matrix["median"].astype(int),
                   super_matrix["mean"].astype(int),
                   s=1, label="Consensus")
        ax.scatter(wrong_matrix["median"].astype(int),
                   wrong_matrix["mean"].astype(int),
                   c="red", s=1, label="Erroneous")
        ax.set_xlabel("Median Gene Length (log)")
        ax.set_ylabel("Mean Gene Length (log)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()
        ax.grid()
        save_or_show(fig, output_dir, "median_vs_mean_length.png")

        # 11. Median length vs gene length
        fig, ax = plt.subplots()
        ax.scatter(super_matrix["median"].astype(int),
                   super_matrix["gene_length"].astype(int),
                   s=1, label="Consensus")
        ax.scatter(wrong_matrix["median"].astype(int),
                   wrong_matrix["gene_length"].astype(int),
                   c="red", s=1, label="Erroneous")
        ax.set_xlabel("Median Gene Length (log)")
        ax.set_ylabel("Gene Length (log)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()
        ax.grid()
        save_or_show(fig, output_dir, "median_vs_gene_length.png")

        # 12. HOG size vs mismatch count
        fig, ax = plt.subplots()
        ax.scatter(wrong_matrix["hogSize"].astype(int),
                   wrong_matrix["mismatchCount"].astype(int),
                   c="blue", s=2)
        ax.set_title("HOG Size vs HMM Disagreement Count")
        ax.set_xlabel("HOG Gene Count")
        ax.set_ylabel("HOG HMM Disagree Gene Count")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid()
        save_or_show(fig, output_dir, "hogsize_vs_mismatch.png")

        # 13. Z-score histogram: all genes
        z_scores = super_matrix["z_scores"]
        z_filtered = np.where(np.isfinite(z_scores), z_scores, 0)
        fig, ax = plt.subplots()
        ax.hist(z_filtered, density=True, bins=100)
        ax.set_title("Z-Score of Gene Length in Parent HOG")
        ax.set_xlabel("Z-score")
        ax.set_ylabel("Frequency")
        ax.grid()
        save_or_show(fig, output_dir, "zscore_histogram_all.png")

        # 14. Z-score histogram: erroneous hits only
        fig, ax = plt.subplots()
        ax.hist(wrong_matrix["z_scores"], density=True, bins=100, color="red")
        ax.set_title("Z-Score of Gene Length — Erroneous Hits")
        ax.set_xlabel("Z-score")
        ax.set_ylabel("Frequency")
        ax.grid()
        save_or_show(fig, output_dir, "zscore_histogram_erroneous.png")

    # 15. Histogram: disagree ratio distribution
    if disagree_df is not None:
        fig, ax = plt.subplots()
        ax.hist(disagree_df["ratio"], bins=100)
        ax.set_title("Ratio of Disagreement Between HMM and HOG")
        ax.set_xlabel("Ratio")
        ax.set_ylabel("Frequency")
        ax.grid()
        save_or_show(fig, output_dir, "disagree_ratio_histogram.png")

    # 18 & 19. IQR elbow plots
    if confusion_df is not None:
        cutoffs = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
        lengths = []

        def calculate_quartiles_and_filter(group, multiplier):
            Q1 = group["Score"].quantile(0.25)
            Q3 = group["Score"].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - multiplier * IQR
            return group[group["Score"] >= lower_bound]

        for cutoff in cutoffs:
            filtered = (
                confusion_df
                .groupby("HmmHogHit", group_keys=False)
                .apply(lambda x: calculate_quartiles_and_filter(x, cutoff))
                .reset_index(drop=True)
            )
            lengths.append(len(filtered))

        big = len(confusion_df)

        # 18. Genes remaining
        fig, ax = plt.subplots()
        ax.plot(cutoffs, lengths, marker="o")
        ax.set_xlabel("IQR Multiplier")
        ax.set_ylabel("Number of Genes Remaining")
        ax.set_title("Genes Remaining by IQR Multiplier")
        ax.axvline(x=iqr_multiplier, color="red", linestyle="--",
                   label=f"Current: {iqr_multiplier}")
        ax.legend()
        ax.grid()
        save_or_show(fig, output_dir, "iqr_elbow_remaining.png")

        # 19. Genes removed
        outs = [big - x for x in lengths]
        fig, ax = plt.subplots()
        ax.plot(cutoffs, outs, marker="o")
        ax.set_xlabel("IQR Multiplier")
        ax.set_ylabel("Number of Genes Removed")
        ax.set_title("Genes Removed by IQR Multiplier")
        ax.axvline(x=iqr_multiplier, color="red", linestyle="--",
                   label=f"Current: {iqr_multiplier}")
        ax.legend()
        ax.grid()
        save_or_show(fig, output_dir, "iqr_elbow_removed.png")


# ---------------------------------------------------------------------------
# Network Plots (from confusionMatrxPlots.ipynb)
# ---------------------------------------------------------------------------

def plot_network(data, output_dir):
    """Generate HOG disagreement network plots."""
    print("Generating network plots...")

    wrong_matrix = data.get("wrong_matrix")
    if wrong_matrix is None:
        print("Warning: wrongMatrix not found. Skipping network plots.")
        return

    # 16. Full disagreement directed graph
    G = nx.DiGraph()
    to_graph = wrong_matrix[["Hog", "HmmHogHit"]].drop_duplicates()
    G.add_edges_from(to_graph[["Hog", "HmmHogHit"]].to_records(index=False))

    fig, ax = plt.subplots(figsize=(12, 12))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, font_size=5, node_size=100, ax=ax)
    ax.set_title("Directed Graph of Incorrect HMM Assignment")
    save_or_show(fig, output_dir, "disagreement_network.png")

    # Build filtered wrong matrix (HOGs that appear as both source and target)
    hog_set = set(wrong_matrix["Hog"])
    hmm_set = set(wrong_matrix["HmmHogHit"])
    intersect = hog_set.intersection(hmm_set)
    filtered = pd.concat([
        wrong_matrix[wrong_matrix["Hog"].isin(intersect)],
        wrong_matrix[wrong_matrix["HmmHogHit"].isin(intersect)]
    ])

    # 17. Cyclic disagreements only
    G2 = nx.DiGraph()
    to_graph2 = filtered[["Hog", "HmmHogHit"]].drop_duplicates()
    G2.add_nodes_from(to_graph2["Hog"].unique())
    G2.add_edges_from(to_graph2[["Hog", "HmmHogHit"]].to_records(index=False))

    cycles = list(nx.simple_cycles(G2))
    cyclic_graph = nx.DiGraph()
    for cycle in cycles:
        cyclic_graph.add_nodes_from(cycle)
        cyclic_graph.add_edges_from(zip(cycle, cycle[1:] + [cycle[0]]))

    fig, ax = plt.subplots(figsize=(12, 12))
    pos = nx.spring_layout(cyclic_graph)
    nx.draw(cyclic_graph, pos, with_labels=True, font_size=5,
            node_size=100, node_color="red", ax=ax)
    ax.set_title("Cyclic Disagreements in HMM Assignment")
    save_or_show(fig, output_dir, "disagreement_network_cyclic.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = get_args()

    run_dir = f"HmmRun{args.iteration}"
    if not os.path.exists(run_dir):
        print(f"Error: Run directory '{run_dir}' not found.")
        sys.exit(1)

    db_path = os.path.join(
        run_dir,
        f"{args.project_name}HmmDatabase",
        f"{args.project_name}fam1"
    )

    output_dir = args.output_dir
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print(f"--- FamBuilder plots ---")
    print(f"Project  : {args.project_name}")
    print(f"Iteration: {args.iteration}")
    print(f"Plots    : {args.plot}")
    print(f"Output   : {output_dir or 'interactive'}")

    # HOG stats plots
    if args.plot in ("all", "hog_stats"):
        hog_stats = load_hog_stats(run_dir)
        if hog_stats is not None:
            # Detect grass filter by checking if BOP/PACMAD species columns exist
            pacmad = ["Svi", "Pvi", "Pha", "Msi", "Sbi", "Zma", "Ecu", "Pau"]
            grass_filter = all(s in hog_stats.columns for s in pacmad)
            plot_hog_stats(hog_stats, output_dir, grass_filter=grass_filter)

    # Confusion matrix plots
    if args.plot in ("all", "confusion"):
        data = load_confusion_data(run_dir, db_path)
        plot_confusion(data, output_dir, args.iqr_multiplier)

    # Network plots
    if args.plot in ("all", "network"):
        data = load_confusion_data(run_dir, db_path)
        plot_network(data, output_dir)

    print("Done.")


if __name__ == "__main__":
    main()