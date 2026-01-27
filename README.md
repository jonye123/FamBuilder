# GrassFam HMM Building Pipeline

This repository contains scripts for building and refining Hidden Markov Models (HMMs) for gene family classification in grasses, particularly focused on PACMAD and BOP clades.

## Overview

The pipeline takes orthogroups from OrthoFinder output, creates multiple sequence alignments, builds HMMs, and iteratively refines these models through a series of filtering and testing steps. The workflow is designed to create high-quality HMM profiles for accurate gene family classification.

## Requirements

- Python 3.9+ with the following packages:
  - BioPython
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - jupyter
- HMMER3 suite (hmmbuild, hmmscan)
- MAFFT aligner
- GNU Parallel
- micromamba or conda environment manager

## Pipeline Structure

```
hmmBuilding/
├── main.sh                    # Main execution script
├── orthofinder_to_fasta.py    # Convert OrthoFinder output to FASTA
├── hogFilter.ipynb            # Filter HOGs based on taxonomic distribution
├── confusionMatrxPlots.ipynb  # Analyze and visualize results
└── bashScripts/
    ├── run_all.sh             # Wrapper for all bash scripts
    ├── maffter.sh             # Run MAFFT for multiple sequence alignment
    ├── hmmbuilder.sh          # Build HMMs from alignments
    ├── confusionMatrixMaker.sh # Generate confusion matrix for evaluation
    └── ...
```

## Usage

### 1. Input Preparation

Start with OrthoFinder output, specifically the `N0.tsv` file which contains hierarchical orthogroups. Include a path to the directory of fastas for each species, titled as {3letterspeciescode}.fasta

### 2. Running the Pipeline

Execute the main script to run the full pipeline:

```bash main.sh {number_threads}
```

This will:
1. Create a new directory for the current iteration
2. Convert orthogroups to FASTA files
3. Filter HOGs based on taxonomic representation
4. Generate multiple sequence alignments
5. Build HMMs
6. Test and evaluate the models
7. Iterate if necessary to improve results

## Key Parameters

In `hogFilter.ipynb`, you can adjust the function hog_filter:
- Minimum species count threshold (default: 8)
- PACMAD and BOP representation requirements

In `hmmbuilder.sh`, you can configure:
- CPU usage for hmmbuild
- Output directory structure

## Output Files

The pipeline produces:
- Filtered HOG FASTA files in `./hogDirectoryFiltered/`
- Multiple sequence alignments in `./MSA/`
- HMM profiles in `./grassfamHmmLib/`
- Evaluation metrics and confusion matrices in `./grassHmmDatabase/grassfam1/`

## Iterative Refinement

The `main.sh` script automatically manages iterations, creating a new directory (`HmmRun{iteration}`) for each run. The process continues until convergence criteria are met.

## Troubleshooting

- Ensure all input FASTA files are in the correct format
- Check permissions for output directories
- Verify that all required software is in your PATH
- For memory issues, adjust the parallelization in bash scripts