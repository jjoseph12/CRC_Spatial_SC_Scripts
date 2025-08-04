For these scripts, adjust the HPC settings/environment for personal use. 

run_stardist.py: Runs StarDist2D nucleus segmentation on a single high-res H&E image. Automatic Parameters based on image size and saves segmentation masks and overlay PNGs.
stardist_batch.sh: SLURM batch job script to run run_stardist.py on all images in a directory. Skips already-processed images and manages memory/resource settings for HPC environments.

scrna_integration_Pelka_Moorman_10x.ipynb: Performs integration of 10x Genomics scRNA-seq data with the Pelka and Moorman single-cell atlas using Harmony, followed by label transfer.

single_cell_analysis.ipynb: Full analysis pipeline for a single-cell dataset, including preprocessing (HVG selection, PCA, Leiden), embedding (UMAP), and marker gene identification (rank_genes_groups).

old_visium_single_cell_integration.ipynb: Older integration workflow between Visium spatial data and single-cell references, older outdated methods or older label mappings.

automated_single_cell.ipynb: Automates the full single-cell label transfer pipeline using scanpy and ingest, for multiple samples. Uses Harmony for batch correction.

label_single_cell.ipynb: Focused notebooks that perform only the labeling portion takes a reference and a query dataset and applies predicted cell types. 
