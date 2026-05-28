# FamBuilder

Iterative HMM-based gene family classification pipeline. Generalized from Noor Mashiana's GrassFam.

## Requirements

The following tools must be available in your PATH:

- HMMER v3.4+ (`hmmbuild`, `hmmscan`, `hmmpress`)
- MAFFT
- GNU Parallel

## Installation

```bash
pip install fambuilder
```

Or for development:

```bash
git clone https://github.com/jonye123/FamBuilder
cd FamBuilder
pip install -e .
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
fambuilder --config my_config.yaml
```

4. To resume after an interruption:

```bash
fambuilder --config my_config.yaml --resume
```

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