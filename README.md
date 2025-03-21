# Under active development

# LIANA-NicheNet Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.3-brightgreen.svg)](https://www.nextflow.io/)
[![GitHub Actions Status](https://github.com/vishaloza/ccc-pipeline/workflows/tests/badge.svg)](https://github.com/vishaloza/ccc-pipeline/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Nextflow pipeline for integrated cell-cell communication analysis using LIANA (Python) and NicheNet (R).

## Quick Start

You can run this pipeline directly from GitHub with Nextflow:

```bash
# Run with test data (requires Nextflow and Conda)
nextflow run yourusername/liana-nichenet-pipeline -r main \
    --input your_data.h5ad \
    --sender_celltype TypeA \
    --receiver_celltype TypeB \
    -profile conda

# Run on an HPC with SLURM
nextflow run yourusername/liana-nichenet-pipeline -r main \
    --input your_data.h5ad \
    --sender_celltype TypeA \
    --receiver_celltype TypeB \
    -profile slurm,singularity
```

## Pipeline Overview

This pipeline integrates:

1. **LIANA**: For ligand-receptor interaction identification
2. **NicheNet**: For downstream gene regulation prediction

The workflow:
- Runs LIANA on AnnData object to prioritize ligand-receptor pairs
- Feeds these pairs into NicheNet for downstream analysis
- Creates integrated visualizations connecting cell communication to gene regulation

## Requirements

- Nextflow (≥21.10.3)
- One of the following:
  - Docker or Singularity (recommended for HPC)
  - Conda/Mamba

## Input Requirements

- An AnnData object (`.h5ad` file) with scRNA-seq data
- Cell type annotations should be in `adata.obs['cell_type']`
- Sender and receiver cell type names (must match names in the data)

## Parameters

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `--input` | Path to input `.h5ad` file | Yes | None |
| `--sender_celltype` | Name of sender cell type | Yes | None |
| `--receiver_celltype` | Name of receiver cell type | Yes | None |
| `--outdir` | Output directory | No | `results` |

## Output

The pipeline produces:

```
results/
├── liana/
│   ├── liana_ranked_interactions.csv   # Ranked L-R interactions
│   ├── top_lr_pairs_for_nichenet.csv   # Top pairs for NicheNet input
│   └── liana_dotplot.pdf               # Visualization
├── nichenet/
│   ├── nichenet_results_from_liana.rds  # R object with results
│   ├── ligand_activity_heatmap.pdf      # Activity predictions 
│   └── ligand_target_heatmap.pdf        # Target gene regulation
└── integrated_results/
    ├── liana_nichenet_integrated_results.csv  # Combined results
    ├── liana_nichenet_correlation.pdf         # Method correlation
    └── top_ligands_summary_heatmap.pdf        # Ranking comparison
```

## Execution Profiles

| Profile | Description |
|---------|-------------|
| `conda` | Uses conda environments defined in `./envs/` |
| `docker` | Uses Docker containers |
| `singularity` | Uses Singularity containers (recommended for HPC) |
| `slurm` | Configuration for SLURM scheduler |
| `sge` | Configuration for SGE scheduler |

Combine profiles with comma: `-profile slurm,singularity`

## Citation

If you use this pipeline, please cite:

- LIANA: https://github.com/saezlab/liana-py
- NicheNet: https://github.com/saeyslab/nichenetr
- Nextflow: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
