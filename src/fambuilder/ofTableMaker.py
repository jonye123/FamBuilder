#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse


def confusion_to_ortho(table_in, table_out, column, header, best, extractspecies, evalue_threshold, no_threshold):

    if column:
        hog_table = pd.read_csv(table_in, delimiter="\t")[["Gene", column]]

        print("Extracting species from Gene names")
        hog_table['Species'] = hog_table["Gene"].apply(lambda x: x.split("|")[0])

        gene_groups = hog_table.groupby([column, 'Species'])['Gene'].apply(lambda x: ','.join(x))
        pivot_df = gene_groups.unstack(fill_value='').reset_index()

    else:
        hog_table = pd.read_csv(table_in, delimiter="\t", header=0 if header else None)
        if not header:
            hog_table.columns = ["Gene", "Species", "HmmHogHit", "E-value", "Score"]
        with np.errstate(divide='ignore', invalid='ignore'):
            hog_table["log_E"] = np.log10(hog_table["E-value"].astype(float))

        # allow HmmHogHit to be used as Orthogroup
        hog_table.rename(columns={"HmmHogHit": "Orthogroup"}, inplace=True)

        if not no_threshold:
            hog_table = hog_table[hog_table["log_E"] <= evalue_threshold]

        if extractspecies:
            print("Extracting species from Gene names")
            hog_table['Species'] = hog_table["Gene"].apply(lambda x: x.split("|")[0])

        if best:
            hog_table = hog_table.sort_values(by="E-value", ascending=True)
            best_hits = hog_table.groupby(['Orthogroup', 'Species'], as_index=False).first()
            pivot_df = best_hits.pivot(index='Orthogroup', columns='Species', values='Gene').reset_index()
        else:
            gene_groups = hog_table.groupby(['Orthogroup', 'Species'])['Gene'].apply(lambda x: ','.join(x))
            pivot_df = gene_groups.unstack(fill_value='').reset_index()

    pivot_df.to_csv(table_out, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert hmmscore confusion matrix to orthofinder format.')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('output', help='Path to the output file')
    parser.add_argument('--column', type=str, default='', help='Column name to use for cluster')
    parser.add_argument('--header', action='store_true', help='Flag to indicate if the input file has a header')
    parser.add_argument('--best', action='store_true', help='Flag to indicate keeping just the best gene match for each hmm/hog')
    parser.add_argument('--extractspecies', dest='extractspecies', action='store_true', help='Flag to indicate if Species needs to be extracted from Gene names')
    parser.add_argument('--evalue_threshold', type=float, default=-23, help='E-value threshold for filtering hits (default: -23)')
    parser.add_argument('--no-threshold', action='store_true', help='Flag to indicate no E-value threshold filtering')

    args = parser.parse_args()
    confusion_to_ortho(args.input, args.output, args.column, args.header, args.best, args.extractspecies, args.evalue_threshold, args.no_threshold)
