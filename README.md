# Overview of Scripts and Notebooks

> **Note:** For all scripts, be sure to adjust HPC settings and environment variables according to your personal cluster configuration.

# Conda Environment Setup for ENACT and Integrated Spatial Transcriptomics Pipelines

This directory contains Conda environment YAML files used to recreate specific software environments for the ENACT pipeline and for integrated spatial transcriptomics tools such as SMURF, Bin2Cell, and general downstream analysis workflows.

## Environment Files

| File Name                      | Purpose                                               | Build Information     | Notes                                                  |
|-------------------------------|-------------------------------------------------------|------------------------|--------------------------------------------------------|
| enact_py_env_full.yml         | Environment used for the ENACT pipeline               | Includes build strings | Full reproducibility across systems                   |
| enact_py_env_nobuilds.yml     | Portable version of the ENACT environment             | No build strings       | More flexible for installing on different machines     |
| integrate_env_full.yml        | Environment used for SMURF, BinCell, and analysis    | Includes build strings | Used for all spatial transcriptomics analysis pipelines |
| integrate_env_nobuilds.yml    | Portable version for SMURF, Bin2Cell, and analysis    | No build strings       | Easier to set up on different systems                 |

---

## Environment Use Cases

- `enact_py_env`: Dedicated to the ENACT pipeline only.
- `integrate_env`: Shared environment used across:
  - SMURF segmentation pipeline
  - Bin2Cell reconstruction
  - General downstream analysis (e.g., Scanpy, visualization)

---

## How to Create a Conda Environment from a YAML File

### Option 1: Use the name defined in the YAML

This will create an environment using the name specified inside the YAML file:

```bash
conda env create -f integrate_env_full.yml

### Option 2: You can replace the default name with one of your choosing by using the -n flag:
conda env create -f integrate_env_full.yml -n smurf_env
conda env create -f enact_py_env_nobuilds.yml -n enact_custom


## How to Create a Conda Environment from a YAML File

conda activate enact_py_env

or if custom name: conda activate enact_custom


## Models
- `/gpfs/commons/groups/innovation/jjoseph/python_2D_versatile_he.zip`
-   Contains the trained StarDist 2D versatile model for H&E images.


---

## Helpful Scripts
- **Launch Jupyter Notebook**: `Jupyter.sh`  
  Starts a Jupyter Notebook session on an HPC cluster.
- **Interactive SMURF Job:** `interactive_smurf.sh`
  - Launches GPU Node for interactive, right now its set for SMURF environment but can be adjusted for differnt purposes. 
- 


## Python Scripts

### `run_stardist.py`
Applies StarDist2D nucleus segmentation to a single high-resolution H&E image.

- Automatically selects parameters based on image dimensions  
- Outputs segmentation masks and overlay PNGs

### `stardist_batch.sh`
SLURM batch script to process a directory of images using `run_stardist.py`.

- Skips already processed images  
- Optimized for memory and resource handling in HPC environments

### ENACT 

### `run_enact_job.sh`
Launches ENACT segmentation and analysis workflow for a given input image or batch. Individual one with given config file

# Now for Processing Multiple at a time: 

## generate_configs.py: 
- Generates configuration files used for ENACT

### `run_generate_config_batch.sh`
Bash script to run `generate_configs.py` across multiple samples.

### `run_individual_sample_batch.sh`
Executes ENACT pipeline for a single sample using batch SLURM submission, you can specify which one, this can be used with submit_all_enact_jobs.sh to run them all. 

### `submit_all_enact_jobs.sh`
Master script to submit multiple ENACT jobs via SLURM.


### Bin2cell 

'run_bin2cell_full.py' = script that runs bin2cell with given data input for all images 
' run_bin2cell_full.py' = bash script to run bin2cell 

## Jupyter Notebooks

Analyis
### 'Qc_Analysis': 
- General analysis script (details not provided in description).

### `analysis-4.ipynb`
Similar analysis additional spatial transcriptomics or QC analysis for specific samples.


## Smurf (there has been some changes to code will talk about) 
### `07_18_p2_smurf-4.ipynb`
SMURF-based spatial reconstruction and analysis notebook for the P2 sample. (Worked all the way through but needs some modifcations)

### `07_25_p5_smurf-2.ipynb`
SMURF-based spatial reconstruction and analysis notebook for the P5 sample. (No changes to code but does get stuck at some steps, more to this later)

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
