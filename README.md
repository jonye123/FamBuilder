# FamBuilder

Iterative HMM-based gene family classification pipeline. Generalized from Noor Mashiana's GrassFam.

## Requirements

The following tools must be available in your PATH:

- HMMER v3.4+ (`hmmbuild`, `hmmscan`, `hmmpress`)
- MAFFT
- GNU Parallel
- seqkit (only required for `scripts/splitAndScore.sh`)

## Installation
For development:

```bash
git clone https://github.com/jonye123/FamBuilder
cd FamBuilder
pip install -e .
```

Once published:
```bash
pip install fambuilder
```

## Usage

1. Create a working directory and place your `N0.tsv` file inside it:

```bash
mkdir my_project
cd my_project
```

2. Copy and edit the config template:

```bash
cp config.yaml my_config.yaml
# edit my_config.yaml with your paths and settings
```

3. Run the pipeline:

```bash
fambuilder run --config my_config.yaml
```

4. To resume after an interruption:

```bash
fambuilder run --config my_config.yaml --resume
```

## Scoring

Once a pipeline run has finished (`HMMFinal/` exists), you can score additional peptide FASTAs — e.g. a new species' proteome — against the final HMM library:

```bash
fambuilder score --config my_config.yaml --source_fasta new_species.fasta
```

This locates the library at `HMMFinal/<project_name>` using `project_name` from the config. To score against a library elsewhere, pass it directly instead:

```bash
fambuilder score --hmm_lib path/to/hmm_lib --source_fasta new_species.fasta --jobs 8
```

Optional flags: `--jobs` (parallel splitting/scanning jobs, default 8), `--split_dir`, `--output_dir`, `--confusion_dir` (default to `./split`, `./hmmerOutput`, `./catThis`).

### Converting scoring output to an OrthoFinder table

`fambuilder score` writes one confusion file per species chunk under `./catThis/` (tab-separated, no header: `Gene`, `Species`, `HmmHogHit`, `E-value`, `Score`). The plotting commands below expect an OrthoFinder-style table instead (an `Orthogroup` column plus one column per species). `ortho-table` converts between the two:

```bash
cat catThis/* > all_confusion.tsv
fambuilder ortho-table all_confusion.tsv Orthogroups.tsv
```

By default it drops any hit with `log10(E-value) > -23`; pass `--no-threshold` to keep everything, or `--evalue_threshold` to set your own cutoff. By default all hits for a given `(Orthogroup, Species)` pair are comma-joined into one cell; pass `--best` to keep only the single lowest-E-value gene instead. Pass `--header` if your input file already has a header row (e.g. `confusionMatrixMaker.sh`'s `confusionMatrix` output), and `--extractspecies` if species should be parsed from the `SpeciesCode|GeneID`-style `Gene` column rather than read from a `Species` column. Run `fambuilder ortho-table --help` for the full flag list.

## Plotting

Three standalone plotting commands work directly on an OrthoFinder table (`Orthogroups.tsv` or a HOG table like `N0.tsv`) — no FamBuilder run required. They require the `plots` extra (`pip install fambuilder[plots]`). Each forwards its arguments to its own script, so run any of them with `--help` for the full flag list.

**`figures`** — genes-per-family and species-per-family histograms:

```bash
fambuilder figures Orthogroups.tsv --clade my_clade --min_species 4
```

Writes `<clade>_genes_per_family_hist.png` and `<clade>_species_per_family_hist.png` (`output_*` if `--clade` is omitted).

**`heatmap`** — species-presence heatmap comparing the two root clades of a Newick tree:

```bash
fambuilder heatmap Orthogroups.tsv tree.nwk --clade1 GroupA --clade2 GroupB
```

Writes `<clade1>_vs_<clade2>_heatmap.png` (or the path given via `-o/--output`).

**`occurrence`** — OrthoVenn3-style presence/absence pattern plot:

```bash
fambuilder occurrence Orthogroups.tsv -n 40 --min_species 2
```

Writes `presence_absence_patterns.png` by default (`-o/--output` to change it); `--order` takes a file of species names (one per line) to fix column order, and `--tsv` optionally dumps the full pattern table.

## Config Options

| Option | Required | Default | Description |
|---|---|---|---|
| `project_name` | Yes | — | Name used for output directories |
| `input` | Yes | — | Path to OrthoFinder `N0.tsv` |
| `fasta_dir` | Yes | — | Directory of per-species proteome FASTAs |
| `min_species` | No | 4 | Minimum species per HOG |
| `grass_filter` | No | false | Use BOP/PACMAD grass-specific filter |
| `jobs` | No | 8 | Parallel jobs for MAFFT/hmmbuild/hmmscan |
| `conda_env` | No | "" | Conda environment to activate |
| `iqr_multiplier` | No | 3.5 | IQR multiplier for score filtering |
| `max_iterations` | No | 15 | Maximum iterations before stopping |
| `convergence_threshold` | No | 0.000001 | Stop when error change falls below this |
| `cleanup` | No | true | Delete intermediate directories on completion |

## Input Format

- `N0.tsv`: OrthoFinder hierarchical orthogroups file
- Proteome FASTAs: one file per species, sequence IDs in `SpeciesCode|GeneID` format

## Output

Final HMM library is written to `HMMFinal/` in your working directory.

## Credits

Built on Noor Mashiana's GrassFam. Generalized by Jon Ye.
