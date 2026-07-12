"""
FamBuilder: Iterative HMM-based gene family classification pipeline.
Generalized from Noor Mashiana's GrassFam.

Usage:
    fambuilder run --config config.yaml
    fambuilder run --config config.yaml --resume
    fambuilder score --config config.yaml --source_fasta new_species.fasta

Modules:
    main                  - Pipeline orchestrator and entry point
    hogFilter             - HOG filtering and FASTA preparation
    confusionMatrixPlots  - HMM confusion matrix analysis
    orthofinder_to_fasta  - OrthoFinder output parser
"""

__version__ = "0.1.0"
__author__ = "Jon Ye"
__license__ = "MIT"