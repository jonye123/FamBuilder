import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse


def count_genes(cell):
    """Count genes in a cell (comma-separated list like 'Ami|123, Ami|456')."""
    if pd.isna(cell) or cell == "":
        return 0
    genes = [g.strip() for g in str(cell).split(",") if g.strip()]
    return len(genes)


def detect_format(columns):
    if "Orthogroup" in columns:
        return "orthogroup"
    elif "HOG" in columns:
        return "hog"
    else:
        raise ValueError("Unknown format: cannot detect Orthofinder table type.")


def get_species_columns(df, fmt):
    if fmt == "orthogroup":
        return [col for col in df.columns if col != "Orthogroup"]
    elif fmt == "hog":
        meta_cols = ["HOG", "OG", "Gene Tree Parent Clade"]
        return [col for col in df.columns if col not in meta_cols]


def process(df):
    fmt = detect_format(df.columns)
    species_cols = get_species_columns(df, fmt)

    genes_per_family = []
    species_per_family = []

    for _, row in df.iterrows():
        gene_counts = [count_genes(row[col]) for col in species_cols]

        total_genes = sum(gene_counts)
        species_count = sum(1 for g in gene_counts if g > 0)

        genes_per_family.append(total_genes)
        species_per_family.append(species_count)

    return genes_per_family, species_per_family


def plot_histogram(data, title, xlabel, output, bins=None, clade=None, xlim=None):
    plt.figure()

    if bins is None:
        bin_size = 1
        bins = np.arange(0, max(data) + bin_size, bin_size)
        plt.hist(data, bins=bins)
    else:
        plt.hist(data, bins=bins)

    if clade:
        title = f"{title} ({clade})"

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")

    if xlim is not None:
        plt.xlim(0, xlim)

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot OrthoFinder histograms")
    parser.add_argument("input_file", help="Orthofinder table (TSV)")
    parser.add_argument("--clade", help="Clade name for labeling", default=None)
    parser.add_argument("--xlim", type=int, default=None,
                        help="Upper x-axis limit for genes-per-family histogram")
    parser.add_argument("--min_species", type=int, default=None,
                        help="Filter genes-per-family histogram to families with >= X species")

    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep="\t", dtype=str)

    genes_per_family, species_per_family = process(df)

    prefix = args.clade if args.clade else "output"

    genes_file = f"{prefix}_genes_per_family_hist.png"
    species_file = f"{prefix}_species_per_family_hist.png"

    # Apply species filter if requested
    if args.min_species is not None:
        filtered_genes = [g for g, s in zip(genes_per_family, species_per_family) if s >= args.min_species]
        genes_title = f"Genes per Family (>= {args.min_species} species)"
    else:
        filtered_genes = genes_per_family
        genes_title = "Genes per Family"

    # Species bins (0.5-centered)
    max_species = max(species_per_family)
    bins = [x + 0.5 for x in range(0, max_species + 1)]

    # Plot genes per family
    plot_histogram(
        filtered_genes,
        genes_title,
        "Number of Genes",
        genes_file,
        clade=args.clade,
        xlim=args.xlim
    )

    # Plot species per family
    plot_histogram(
        species_per_family,
        "Species per Family",
        "Number of Species",
        species_file,
        bins=bins,
        clade=args.clade
    )

    print("Done. Output files:")
    print(f" - {genes_file}")
    print(f" - {species_file}")


if __name__ == "__main__":
    main()
