#!/usr/bin/env python3
"""
Generate an OrthoVenn3-style presence/absence pattern plot from OrthoFinder output.

Each row represents a unique combination of species presence/absence across clusters.
Columns are species. Cells show rounded green bars (present) or faded (absent).
A horizontal bar on the right shows how many clusters match that pattern.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib

matplotlib.rcParams['font.family'] = 'sans-serif'


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
            "Expected 'Orthogroup' or 'HOG' column."
        )


def get_species_columns(df, fmt):
    if fmt == "orthogroup":
        return [col for col in df.columns if col != "Orthogroup"]
    elif fmt == "hog":
        meta_cols = {"HOG", "OG", "Gene Tree Parent Clade"}
        return [col for col in df.columns if col not in meta_cols]


def read_order_file(order_path):
    """Read a newline-separated file of species names defining column order."""
    with open(order_path, "r") as fh:
        names = [line.strip() for line in fh if line.strip()]
    return names


# ---------------------------------------------------------------------------
# Build presence/absence patterns
# ---------------------------------------------------------------------------

def build_pattern_counts(df, fmt, species_cols):
    """
    For each cluster, build a binary presence/absence tuple across species.
    Return a DataFrame with columns = species + 'count', one row per unique
    pattern, sorted by count descending.
    """
    patterns = []
    for _, row in df.iterrows():
        pattern = tuple(1 if count_genes(row[sp]) > 0 else 0 for sp in species_cols)
        patterns.append(pattern)

    pattern_series = pd.Series(patterns)
    counts = pattern_series.value_counts().reset_index()
    counts.columns = ["pattern", "count"]

    rows = []
    for _, r in counts.iterrows():
        rec = {sp: val for sp, val in zip(species_cols, r["pattern"])}
        rec["count"] = r["count"]
        rows.append(rec)

    result = pd.DataFrame(rows)
    result = result.sort_values("count", ascending=False).reset_index(drop=True)
    return result


# ---------------------------------------------------------------------------
# OrthoVenn3-style plot
# ---------------------------------------------------------------------------

def draw_rounded_bar(ax, x, y, width, height, color, alpha=1.0):
    """Draw a rounded rectangle (pill shape) on the given axes."""
    rounding = min(height, width) * 0.4
    fancy = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle=f"round,pad=0,rounding_size={rounding}",
        facecolor=color,
        edgecolor='none',
        alpha=alpha,
        zorder=2,
    )
    ax.add_patch(fancy)


def plot_orthovenn3_style(
    pattern_df,
    species_cols,
    n_rows,
    output_file,
    title=None,
):
    """
    Draw an OrthoVenn3-style presence/absence pattern plot.

    Left panel: pill-shaped bars indicating species presence/absence per pattern row.
    Right panel: horizontal bar chart showing the cluster count per pattern.
    """
    plot_df = pattern_df.head(n_rows).copy().reset_index(drop=True)
    n_patterns = len(plot_df)
    n_species = len(species_cols)

    # --- Figure layout -------------------------------------------------------
    row_height = 0.38
    fig_height = max(6, n_patterns * row_height + 2.0)

    presence_width = n_species * 0.9 + 0.5
    cluster_bar_width = 3.0
    gap = 0.9

    total_width = presence_width + gap + cluster_bar_width + 1.0

    fig = plt.figure(figsize=(total_width, fig_height))
    fig.patch.set_facecolor('#E8F5E9')

    # Manual axes positioning
    left_margin = 0.06
    right_margin = 0.03
    top_margin = 0.06
    bottom_margin = 0.03

    avail_w = 1.0 - left_margin - right_margin
    avail_h = 1.0 - top_margin - bottom_margin

    total_content = presence_width + gap + cluster_bar_width
    frac_presence = presence_width / total_content
    frac_gap = gap / total_content
    frac_cluster = cluster_bar_width / total_content

    ax_presence = fig.add_axes([
        left_margin,
        bottom_margin,
        avail_w * frac_presence,
        avail_h
    ])

    ax_cluster = fig.add_axes([
        left_margin + avail_w * (frac_presence + frac_gap),
        bottom_margin,
        avail_w * frac_cluster,
        avail_h
    ])

    

    # --- Presence/absence matrix (left panel) --------------------------------
    ax_presence.set_facecolor('#E8F5E9')

    bar_width = 0.75
    bar_height = 0.45

    present_color = "#2E8B2E"
    absent_color = "#B5D8B5"

    for pi in range(n_patterns):
        for si in range(n_species):
            val = plot_df.iloc[pi][species_cols[si]]
            if val == 1:
                color = present_color
                alpha = 1.0
            else:
                color = absent_color
                alpha = 0.5

            draw_rounded_bar(
                ax_presence,
                x=si,
                y=pi,
                width=bar_width,
                height=bar_height,
                color=color,
                alpha=alpha,
            )

    ax_presence.set_xlim(-0.6, n_species - 0.4)
    ax_presence.set_ylim(n_patterns - 0.5, -0.5)

    ax_presence.set_xticks(range(n_species))
    ax_presence.set_xticklabels(species_cols, fontsize=24, rotation=45, ha='left', rotation_mode='anchor')
    ax_presence.xaxis.set_ticks_position('top')
    ax_presence.xaxis.set_label_position('top')

    ax_presence.set_yticks([])
    for spine in ax_presence.spines.values():
        spine.set_visible(False)
    ax_presence.tick_params(axis='x', length=0)

    # --- Cluster count bars (right panel) ------------------------------------
    ax_cluster.set_facecolor('#E8F5E9')

    counts = plot_df["count"].values
    max_count = max(counts) if len(counts) > 0 else 1

    for pi in range(n_patterns):
        c = counts[pi]
        intensity = c / max_count
        bar_color = plt.cm.Blues(0.3 + 0.7 * intensity)

        ax_cluster.barh(
            pi, c,
            height=0.6,
            color=bar_color,
            edgecolor='none',
            zorder=2,
        )

        # Count label to the left of bar
        ax_cluster.text(
            -max_count * 0.02, pi,
            str(c),
            ha='right', va='center',
            fontsize=18, fontweight='bold',
            color='#333333',
        )

    ax_cluster.set_xlim(0, max_count * 1.15)
    ax_cluster.set_ylim(n_patterns - 0.5, -0.5)
    ax_cluster.set_yticks([])
    ax_cluster.set_xticks([])

    for spine in ax_cluster.spines.values():
        spine.set_visible(False)

    ax_cluster.set_title("Cluster count", fontsize=24, fontweight='bold', pad=8)

    # --- Optional overall title ----------------------------------------------
    if title:
        fig.suptitle(title, fontsize=13, fontweight='bold', y=0.99)

    plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Plot saved to {output_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate an OrthoVenn3-style presence/absence pattern plot "
            "from OrthoFinder output (Orthogroups.tsv or HOG table)."
        )
    )
    parser.add_argument("input_file", help="OrthoFinder table (TSV)")
    parser.add_argument(
        "-n", "--nrows",
        type=int,
        default=40,
        help="Number of top patterns to display (default: 40)",
    )
    parser.add_argument(
        "-o", "--output",
        default="presence_absence_patterns.png",
        help="Output image filename (default: presence_absence_patterns.png)",
    )
    parser.add_argument(
        "--min_species",
        type=int,
        default=None,
        help="Minimum number of species present in a pattern to be included (default: no filter)",
    )
    parser.add_argument(
        "--max_species",
        type=int,
        default=None,
        help="Maximum number of species present in a pattern to be included (default: no filter)",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom plot title",
    )
    parser.add_argument(
        "--tsv",
        default=None,
        help="Optional: also save the full pattern table to this TSV path",
    )
    parser.add_argument(
        "--order",
        default=None,
        help=(
            "Optional: path to a text file with one species name per line, "
            "defining the column (species) order in the plot"
        ),
    )

    args = parser.parse_args()

    # --- Read input ----------------------------------------------------------
    df = pd.read_csv(args.input_file, sep="\t", dtype=str)
    fmt = detect_format(df.columns)
    species_cols = get_species_columns(df, fmt)

    print(f"Detected format: {fmt}")
    print(f"Species found in table ({len(species_cols)}): {species_cols}")

    # --- Apply custom species order if provided ------------------------------
    if args.order:
        ordered_names = read_order_file(args.order)
        print(f"Species order from file ({len(ordered_names)}): {ordered_names}")

        # Validate: every name in the order file must exist in the table
        missing_in_table = [s for s in ordered_names if s not in species_cols]
        if missing_in_table:
            raise ValueError(
                f"The following species from the order file are not in the "
                f"input table: {missing_in_table}\n"
                f"Available species: {species_cols}"
            )

        # Warn about species in the table but not in the order file
        missing_in_order = [s for s in species_cols if s not in ordered_names]
        if missing_in_order:
            print(
                f"WARNING: The following species are in the table but not in "
                f"the order file and will be excluded: {missing_in_order}"
            )

        # Use only the species listed in the order file, in that order
        species_cols = [s for s in ordered_names if s in species_cols]
        print(f"Final species order ({len(species_cols)}): {species_cols}")

    # --- Build patterns ------------------------------------------------------
    pattern_df = build_pattern_counts(df, fmt, species_cols)

    # --- Apply min_species filter --------------------------------------------
    if args.min_species is not None:
        # Sum the presence (1/0) columns for each pattern row
        pattern_df["_n_species"] = pattern_df[species_cols].sum(axis=1)
        before = len(pattern_df)
        pattern_df = pattern_df[pattern_df["_n_species"] >= args.min_species].copy()
        pattern_df = pattern_df.drop(columns=["_n_species"]).reset_index(drop=True)
        after = len(pattern_df)
        print()
        print()
        print(f"Number of clusters with >= {args.min_species}:     {pattern_df['count'].sum()}")
        print(f"min_species filter ({args.min_species}): {before} → {after} patterns retained")
        print()
        print()

     # --- Apply max_species filter --------------------------------------------
    if args.max_species is not None:
        # Sum the presence (1/0) columns for each pattern row
        pattern_df["_n_species"] = pattern_df[species_cols].sum(axis=1)
        before = len(pattern_df)
        pattern_df = pattern_df[pattern_df["_n_species"] <= args.max_species].copy()
        pattern_df = pattern_df.drop(columns=["_n_species"]).reset_index(drop=True)
        after = len(pattern_df)
        print()
        print()
        print(f"Number of clusters with <= {args.max_species}:     {pattern_df['count'].sum()}")
        print(f"min_species filter ({args.max_species}): {before} → {after} patterns retained")
        print()
        print()

    total_patterns = len(pattern_df)
    print(f"Total unique presence/absence patterns: {total_patterns}")

    n_display = min(args.nrows, total_patterns)
    print(f"Displaying top {n_display} patterns")

    # --- Optionally save full table ------------------------------------------
    if args.tsv:
        out_df = pattern_df.copy()
        out_df.insert(
            0, "pattern_label",
            [
                "".join(str(int(v)) for v in out_df.iloc[i][species_cols])
                for i in range(len(out_df))
            ],
        )
        out_df.to_csv(args.tsv, sep="\t", index=False)
        print(f"Full pattern table saved to {args.tsv}")

    # --- Plot ----------------------------------------------------------------
    plot_orthovenn3_style(
        pattern_df,
        species_cols,
        n_display,
        args.output,
        title=args.title,
    )


if __name__ == "__main__":
    main()