#!/usr/bin/env python3
"""
Parse OrthoFinder output (Orthogroups or HOGs) and a Newick tree to produce
a heatmap of species-presence counts split by the two root clades.
"""

import sys
import argparse
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------------------------
# Newick parsing
# ---------------------------------------------------------------------------

def strip_branch_lengths(newick):
    """Remove branch lengths from a Newick string."""
    return re.sub(r":[0-9eE.\-+]+", "", newick)


def parse_newick_to_nested(newick):
    """
    Parse a Newick string into a nested list structure.
    Each internal node becomes a list of its children; each leaf is a string.
    Handles branch lengths, internal node labels, and quoted labels.
    """
    newick = strip_branch_lengths(newick.strip())

    stack = []
    current = []
    i = 0
    # Track whether we just closed a parenthesis — the next bare token
    # would be an internal-node label rather than a leaf.
    just_closed = False

    while i < len(newick):
        ch = newick[i]

        if ch == "(":
            stack.append(current)
            current = []
            just_closed = False
            i += 1

        elif ch == ")":
            finished = current
            current = stack.pop()
            current.append(finished)
            just_closed = True
            i += 1

        elif ch == ",":
            just_closed = False
            i += 1

        elif ch == ";":
            i += 1

        elif ch in ("'", '"'):
            # Quoted label
            quote = ch
            i += 1
            label = ""
            while i < len(newick) and newick[i] != quote:
                label += newick[i]
                i += 1
            i += 1  # skip closing quote
            if just_closed:
                # Internal node label — discard
                just_closed = False
            else:
                current.append(label)
                just_closed = False

        else:
            # Bare (unquoted) token
            label = ""
            while i < len(newick) and newick[i] not in ("(", ")", ",", ";", ":", "'", '"'):
                label += newick[i]
                i += 1
            label = label.strip()
            if label:
                if just_closed:
                    # Internal node label — discard
                    just_closed = False
                else:
                    current.append(label)
                    just_closed = False

    if len(current) == 1 and isinstance(current[0], list):
        return current[0]
    return current


def collect_leaves(node):
    """Recursively collect all leaf names from the nested-list tree."""
    if isinstance(node, str):
        return [node]
    leaves = []
    for child in node:
        leaves.extend(collect_leaves(child))
    return leaves


def split_root_clades(newick):
    """
    Split the species in a Newick tree into two groups corresponding to
    the two children of the root.  Returns (clade1_species, clade2_species).
    """
    tree = parse_newick_to_nested(newick)

    if len(tree) < 2:
        raise ValueError(
            "Tree root has fewer than 2 children — cannot split into two clades."
        )

    clade1 = collect_leaves(tree[0])
    clade2 = []
    for part in tree[1:]:
        clade2.extend(collect_leaves(part))

    return clade1, clade2


# ---------------------------------------------------------------------------
# OrthoFinder table parsing
# ---------------------------------------------------------------------------

def count_genes(cell):
    """Count genes in a cell (comma-separated list)."""
    if pd.isna(cell) or str(cell).strip() == "":
        return 0
    genes = [g.strip() for g in str(cell).split(",") if g.strip()]
    return len(genes)


def detect_format(columns):
    if "Orthogroup" in columns:
        return "orthogroup"
    elif "HOG" in columns:
        return "hog"
    else:
        raise ValueError(
            "Unknown format: cannot detect OrthoFinder table type. "
            "Expected 'Orthogroup' or 'HOG' column0."
        )


def get_species_columns(df, fmt):
    if fmt == "orthogroup":
        return [col for col in df.columns if col != "Orthogroup"]
    elif fmt == "hog":
        meta_cols = {"HOG", "OG", "Gene Tree Parent Clade"}
        return [col for col in df.columns if col not in meta_cols]


def get_cluster_id_column(fmt):
    return "Orthogroup" if fmt == "orthogroup" else "HOG"


# ---------------------------------------------------------------------------
# Build presence matrix & heatmap
# ---------------------------------------------------------------------------

def build_presence_counts(df, fmt, species_cols):
    """
    Return a DataFrame with one row per cluster and one column per species,
    containing the gene count (converted to presence/absence = 1/0).
    """
    cluster_col = get_cluster_id_column(fmt)
    records = []
    for _, row in df.iterrows():
        rec = {cluster_col: row[cluster_col]}
        for sp in species_cols:
            rec[sp] = 1 if count_genes(row[sp]) > 0 else 0
        records.append(rec)
    return pd.DataFrame(records)


def make_heatmap(
    presence_df,
    clade1_species,
    clade2_species,
    clade1_name,
    clade2_name,
    cluster_col,
    output_file,
):
    """
    Create and save a heatmap: x-axis = clade2 species count per cluster,
    y-axis = clade1 species count per cluster.
    """
    # Ensure columns exist (fill missing species with 0)
    for sp in clade1_species + clade2_species:
        if sp not in presence_df.columns:
            presence_df[sp] = 0

    presence_df[clade1_name] = presence_df[clade1_species].sum(axis=1)
    presence_df[clade2_name] = presence_df[clade2_species].sum(axis=1)

    max_c1 = len(clade1_species)
    max_c2 = len(clade2_species)

    frequency_table = (
        presence_df
        .groupby([clade1_name, clade2_name])
        .size()
        .unstack(fill_value=0)
        .reindex(
            index=range(max_c1 + 1),
            columns=range(max_c2 + 1),
            fill_value=0,
        )
        .sort_index(axis=0, ascending=False)
        .sort_index(axis=1)
    )

    # Custom colormap: white → red
    cmap = LinearSegmentedColormap.from_list("white_red", ["white", "red"])

    fig_size = max(6, (max(max_c1, max_c2) + 1) * 0.7)
    plt.figure(figsize=(fig_size, fig_size), dpi=300)
    ax = sns.heatmap(
        frequency_table,
        cmap=cmap,
        annot=True,
        fmt="d",
        annot_kws={"color": "black"},
        linewidths=0.5,
        linecolor="black",
        cbar=False,
        square=True,
    )
    ax.set_title(
        f"Heatmap of {clade1_name} vs {clade2_name} Species Counts per Cluster",
        fontsize=10,
        pad=15,
    )
    ax.set_xlabel(f"{clade2_name} Count", fontsize=9)
    ax.set_ylabel(f"{clade1_name} Count", fontsize=9)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Heatmap saved to {output_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse OrthoFinder output (Orthogroups.tsv or N0.tsv / HOG table) "
            "and a Newick tree, then produce a heatmap of species-presence "
            "counts split by the two root clades."
        )
    )
    parser.add_argument(
        "input_file",
        help="OrthoFinder table (TSV): Orthogroups.tsv or an N0.tsv / HOG table.",
    )
    parser.add_argument(
        "tree_file",
        help="Newick tree file. May or may not contain branch lengths.",
    )
    parser.add_argument(
        "--clade1",
        default=None,
        help="Display name for the first root clade (default: Clade1).",
    )
    parser.add_argument(
        "--clade2",
        default=None,
        help="Display name for the second root clade (default: Clade2).",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output PNG filename (default: <clade1>_vs_<clade2>_heatmap.png).",
    )

    args = parser.parse_args()

    # --- Read tree -----------------------------------------------------------
    with open(args.tree_file) as fh:
        newick = fh.read().strip()

    clade1_species, clade2_species = split_root_clades(newick)

    clade1_name = args.clade1 if args.clade1 else "Clade1"
    clade2_name = args.clade2 if args.clade2 else "Clade2"

    print(f"{clade1_name} species ({len(clade1_species)}): {clade1_species}")
    print(f"{clade2_name} species ({len(clade2_species)}): {clade2_species}")

    # --- Read OrthoFinder table ---------------------------------------------
    df = pd.read_csv(args.input_file, sep="\t", dtype=str)
    fmt = detect_format(df.columns)
    print(f"Detected format: {fmt}")

    species_cols = get_species_columns(df, fmt)

    # Warn about species in tree but missing from the table and vice versa
    tree_species = set(clade1_species + clade2_species)
    table_species = set(species_cols)

    missing_from_table = tree_species - table_species
    missing_from_tree = table_species - tree_species

    if missing_from_table:
        print(
            f"WARNING: species in tree but not in table (will be all-zero): "
            f"{sorted(missing_from_table)}"
        )
    if missing_from_tree:
        print(
            f"WARNING: species in table but not in tree (ignored in heatmap): "
            f"{sorted(missing_from_tree)}"
        )

    # --- Build presence matrix -----------------------------------------------
    cluster_col = get_cluster_id_column(fmt)
    presence_df = build_presence_counts(df, fmt, species_cols)

    # --- Plot ----------------------------------------------------------------
    output_file = args.output
    if output_file is None:
        output_file = f"{clade1_name}_vs_{clade2_name}_heatmap.png"

    make_heatmap(
        presence_df,
        clade1_species,
        clade2_species,
        clade1_name,
        clade2_name,
        cluster_col,
        output_file,
    )


if __name__ == "__main__":
    main()
