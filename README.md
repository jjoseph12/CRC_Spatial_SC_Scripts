# Overview of Scripts and Notebooks

> **Note:** For all scripts, be sure to adjust HPC settings and environment variables according to your personal cluster configuration.

## Python Scripts

### `run_stardist.py`
Applies StarDist2D nucleus segmentation to a single high-resolution H&E image.

- Automatically selects parameters based on image dimensions  
- Outputs segmentation masks and overlay PNGs

### `stardist_batch.sh`
SLURM batch script to process a directory of images using `run_stardist.py`.

- Skips already processed images  
- Optimized for memory and resource handling in HPC environments

## Jupyter Notebooks

### `scrna_integration_Pelka_Moorman_10x.ipynb`
Integrates 10x Genomics scRNA-seq data with the Pelka & Moorman single-cell reference atlas using Harmony.

- Includes label transfer for annotation

### `single_cell_analysis.ipynb`
Complete single-cell RNA-seq analysis pipeline:

- Preprocessing: HVG selection, PCA, Leiden clustering  
- Embedding: UMAP  
- Marker gene identification: `rank_genes_groups`

### `old_visium_single_cell_integration.ipynb`
Legacy workflow for integrating Visium spatial transcriptomics data with single-cell references.

- Uses outdated methods and/or label mappings  
- Retained for reproducibility or comparison

### `automated_single_cell.ipynb`
Automates label transfer across multiple samples.

- Batch correction with Harmony  
- Uses `scanpy` and `ingest` for scalable integration

### `label_single_cell.ipynb`
Performs cell type annotation by transferring labels from a reference to a query dataset.

- Focused notebook for labeling only
