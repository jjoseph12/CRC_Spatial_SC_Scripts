# Spatial Transcriptomics Pipeline Documentation

> **Important:** Adjust HPC settings and environment variables according to your personal cluster configuration for all scripts.
> More details and better reoganization of files coming as well as more details regading changes internally to code (for enact and smurf)
> as well as location to models for more convenience. 

## Table of Contents
- [Conda Environment Setup](#conda-environment-setup)
- [Models](#models)
- [Utility Scripts](#utility-scripts)
- [Processing Pipelines](#processing-pipelines)
- [Analysis Notebooks](#analysis-notebooks)

---

## Conda Environment Setup

This directory contains Conda environment YAML files for reproducing software environments used in the ENACT pipeline and integrated spatial transcriptomics tools (SMURF, Bin2Cell, downstream analysis).

### Environment Files

| File Name | Purpose | Build Information | Use Case |
|-----------|---------|-------------------|----------|
| `enact_py_env_full.yml` | ENACT pipeline environment | Includes build strings | Full reproducibility across systems |
| `enact_py_env_nobuilds.yml` | Portable ENACT environment | No build strings | Flexible installation on different machines |
| `integrate_env_full.yml` | SMURF, Bin2Cell, and analysis | Includes build strings | Complete spatial transcriptomics analysis |
| `integrate_env_nobuilds.yml` | Portable integration environment | No build strings | Easy setup on different systems |

### Environment Use Cases

- **`enact_py_env`**: Dedicated to the ENACT pipeline only
- **`integrate_env`**: Shared environment for:
  - SMURF segmentation pipeline
  - Bin2Cell reconstruction
  - General downstream analysis (Scanpy, visualization)

### Creating Environments

#### Option 1: Use default name from YAML
```bash
conda env create -f integrate_env_full.yml
```

#### Option 2: Custom environment name
```bash
conda env create -f integrate_env_full.yml -n smurf_env
conda env create -f enact_py_env_nobuilds.yml -n enact_custom
```

#### Activating Environments
```bash
# Default name
conda activate enact_py_env

# Custom name
conda activate enact_custom
```

---

## Models

### StarDist Model
- **Path**: `/gpfs/commons/groups/innovation/jjoseph/python_2D_versatile_he.zip`
- **Description**: Trained StarDist 2D versatile model for H&E image nucleus segmentation

---

## Utility Scripts

### System Management
- **`Jupyter.sh`**: Launch Jupyter Notebook session on HPC cluster
- **`interactive_smurf.sh`**: Launch GPU node for interactive work (configured for SMURF environment, adjustable)

---

## Processing Pipelines

### StarDist Nucleus Segmentation

StarDist nucleus segmentation adapted from 10x Genomics Visium HD segmentation guide. SMURF notebooks also perform StarDist segmentation in their initial sections.

https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd

#### `run_stardist.py`
Single image nucleus segmentation using StarDist2D.
> **Note**: Automatic parameter selection based on image dimensions, outputs segmentation masks and overlay PNGs

#### `stardist_batch.sh`
SLURM batch script for directory-wide StarDist processing.
> **Note**: Skips already processed images, optimized for HPC memory and resource handling

### ENACT Pipeline

#### Single Sample Processing
- **`run_enact_job.sh`**: Launch ENACT workflow for individual samples with configuration file

#### Batch Processing
- **`generate_configs.py`**: Generate ENACT configuration files
- **`run_generate_config_batch.sh`**: Batch configuration file generation
- **`run_individual_sample_batch.sh`**: Execute ENACT for single sample via SLURM
- **`submit_all_enact_jobs.sh`**: Master script for submitting multiple ENACT jobs

### Bin2Cell Pipeline

- **`run_bin2cell_full.py`**: Complete Bin2Cell processing script for all images
- **`run_bin2cell_full.sh`**: Bash wrapper for Bin2Cell execution

---

## Analysis Notebooks

### Quality Control and General Analysis

#### `Qc_Analysis.ipynb`
General quality control and analysis workflows.

#### `analysis-4.ipynb`
Additional spatial transcriptomics and QC analysis for specific samples.

### SMURF Analysis

> **Note:** SMURF code has undergone recent changes and may require modifications.

#### `07_18_p2_smurf-4.ipynb`
SMURF spatial reconstruction and analysis for P2 sample.
> **Note**: Complete workflow but needs modifications

#### `07_25_p5_smurf-2.ipynb`  
SMURF spatial reconstruction and analysis for P5 sample.
> **Note**: No code changes needed but encounters processing bottlenecks

### Spatial Analysis

#### `Spatial_Clustering.ipynb`
Performs preprocessing, normalization, and integration for spatial data.
> **Note**: Uses Harmony for integration, applies label transfer from single-cell reference to spatial data, outputs annotated datasets for downstream analysis

#### `Spatial_Clustering_code.ipynb`
Refactored core functions from `Spatial_Clustering.ipynb` into reusable components.
> **Note**: Contains batch processing and plotting functions for UMAPs and label assignments, saves results and visualizations to organized output directories

#### `graphing_clusters.ipynb`
Visualizes clustering results and analyzes cluster compositions.
> **Note**: Creates UMAPs, bar plots, dot plots, computes cluster purity metrics using k-nearest-neighbor analysis, generates publication-quality plots

#### `single_cell_analysis.ipynb`
Complete scRNA-seq analysis pipeline.
> **Note**: Includes HVG selection, PCA, Leiden clustering, UMAP generation, and marker gene identification

#### `scrna_integration_Pelka_Moorman_10x.ipynb`
Integration of 10x Genomics scRNA-seq data with Pelka & Moorman reference atlas.
> **Note**: Uses Harmony with label transfer for cell type annotation

#### `automated_single_cell.ipynb`
Automated label transfer across multiple samples.
> **Note**: Uses Harmony batch correction and Scanpy ingest for scalable integration

#### `label_single_cell.ipynb`
Focused cell type annotation workflow.
> **Note**: Transfers labels from reference to query datasets

### Single-Cell Analysis

#### `old_visium_single_cell_integration.ipynb`
Legacy Visium spatial transcriptomics integration workflow.
> **Note**: Uses outdated methods/mappings, retained for reproducibility and comparison

### Legacy Analysis

1. **Set up environment**: Choose appropriate YAML file and create conda environment
2. **Configure HPC settings**: Adjust batch scripts for your cluster
3. **Download models**: Ensure StarDist model is accessible at specified path
4. **Run pipeline**: Start with StarDist segmentation, then proceed to ENACT or SMURF
5. **Analyze results**: Use provided notebooks for downstream analysis

## Support

For issues with specific pipelines or notebooks, refer to the individual script documentation and ensure all dependencies are properly installed in the correct conda environment.
