#!/usr/bin/env python
"""
Full Bin2Cell pipeline with every step in the tutorial notebook,
with the QC visualisations.

There are some changes to the code for example, I have to comment some filtering because some datasets wouldn't 
work because of to low counts. 
How to run: 
-----
    python run_bin2cell_full.py P2CRC [--qc]    # create QC PNGs
    python run_bin2cell_full.py P5CRC           # skip QC

For individual runs in interactive for batch runs

Use with submit_b2c.slurm, you can submit both all together or indivudal in batch too 
"""

import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "50000000000" # 50B pixels otherwise it crashes

import sys, argparse, os, cv2, numpy as np, scanpy as sc, bin2cell as b2c
import matplotlib; matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from pathlib import Path
import scipy
import scipy.sparse
import datetime

debug_logs = []

def log_debug(msg):
    print(msg)
    debug_logs.append(f"{datetime.datetime.now().isoformat()} - {msg}")

#Edit this dictionary to add / move samples
SAMPLES = {
    "P2CRC": dict(
        path002    = Path("/gpfs/commons/groups/innovation/jjoseph/data/P2_CRC/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/jjoseph/data/P2_CRC/inputs/Visium_HD_Human_Colon_Cancer_P2_tissue_image.btf"),
        sr_spatial = Path("/gpfs/commons/groups/innovation/jjoseph/data/P2_CRC/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/P2CRC/"),
        shift      = (0, -5),         
    ),
    "P5CRC": dict(
        path002    = Path("/gpfs/commons/groups/innovation/jjoseph/data/P5_CRC/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/jjoseph/data/P5_CRC/inputs/Visium_HD_Human_Colon_Cancer_P5_tissue_image.btf"),
        sr_spatial = Path("/gpfs/commons/groups/innovation/jjoseph/data/P5_CRC/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/P5CRC/"),
        shift      = (0, -5),
    ),
        "P1_CRC": dict(
        path002    = Path("/gpfs/commons/groups/innovation/jjoseph/data/P1_CRC/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/jjoseph/data/P1_CRC/inputs/Visium_HD_Human_Colon_Cancer_P1_tissue_image.btf"),
        sr_spatial = Path("/gpfs/commons/groups/innovation/jjoseph/data/P1_CRC/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/P1_CRC/"),
        shift      = (0, -5),
    ),
    "P3_NAT": dict(
        path002    = Path("/gpfs/commons/groups/innovation/jjoseph/data/P3_NAT/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/jjoseph/data/P3_NAT/inputs/Visium_HD_Human_Colon_Normal_P3_tissue_image.btf"),
        sr_spatial = Path("/gpfs/commons/groups/innovation/jjoseph/data/P3_NAT/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/P3_NAT/"),
        shift      = (0, -5),
    ),
    "P5_NAT": dict(
        path002    = Path("/gpfs/commons/groups/innovation/jjoseph/data/P5_NAT/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/jjoseph/data/P5_NAT/inputs/Visium_HD_Human_Colon_Normal_P5_tissue_image.btf"),
        sr_spatial = Path("/gpfs/commons/groups/innovation/jjoseph/data/P5_NAT/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/P5_NAT/"),
        shift      = (0, -5),
    ),
        "JGCRC14NT": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC14NT/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/CRC14NT_CRC5T/highres/FC_20240808_JGCRC14NT_HE_whole-Spot000001.tif"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC14NT/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC14NT/"),
        shift      = (0, -5),
    ),
    "JGCRC5T": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC5T/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/CRC14NT_CRC5T/highres/FC_20240808_JGCRC5T_HE_whole-Spot000001+A.tif"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC5T/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC5T/"),
        shift      = (0, -5),
    ),
    "JGCRC27T": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC27T/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/CRC27T_CRC41T/Metafer/FC-20240916-VHD-CRC27T-HE_whole-Spots/FC-20240916-VHD-CRC27T-HE_whole-Spot000001.tif"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC27T/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC27T/"),
        shift      = (0, -5),
    ),
    "JGCRC41T": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC41T/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/CRC27T_CRC41T/Metafer/FC-20240916-VHD-CRC41T-HE_whole-Spots/CRC-20240916-VHD-CRC41T-HE_whole-Spot000001.tif"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC41T/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC41T/"),
        shift      = (0, -5),
    ),
    "JGCRC47TA": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC47TA/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/CRC47T/SL-SF-20241204-VisiumHD-JGCRC-47a-Spot000005.btf"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC47TA/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC47TA/"),
        shift      = (0, -5),
    ),
    "JGCRC47TB": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC47TB/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/CRC47T/SL-SF-20241204-VisiumHD-JGCRC-47b-Spot000003.btf"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC47TB/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC47TB/"),
        shift      = (0, -5),
    ),
    "JGCRC46": dict(
        path002    = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC46/outs/binned_outputs/square_002um/"),
        source_he  = Path("/gpfs/commons/groups/innovation/share/Metafer/202503_SF_ET/SF-ET-20250313-VisiumHD-JGCRC46-Spots/SF-ET-20250313-VisiumHD-JGCRC46-Spot000001.btf"),
        sr_spatial = Path("/gpfs/commons/home/slee/innovation/share/Spaceranger/ultima/JGCRC46/outs/spatial/"),
        outdir     = Path("/gpfs/commons/groups/innovation/jjoseph/bin2cell_results/JGCRC46/"),
        shift      = (0, -5),
    ),
}

#   parameters  (identical to notebook defaults)
MICS_PER_PIXEL     = 0.30
PROB_HE            = 0.10
PROB_GEX           = 0.01
MAX_BIN_DISTANCE   = 4
THR_COUNTS         = 10        

# QC crop windows used in the notebook
MASK1 = dict(r=(2050, 2250), c=(1350, 1550))
MASK2 = dict(r=(2225, 2275), c=(1400, 1450))

def subset(ad, mask):
    return ad[(ad.obs['array_row'].between(*mask["r"])) &
              (ad.obs['array_col'].between(*mask["c"]))].copy()

def run(sample_key: str, do_qc: bool):
    cfg = SAMPLES[sample_key]
    cfg["outdir"].mkdir(parents=True, exist_ok=True)
    os.makedirs(cfg["outdir"]/ "stardist", exist_ok=True)
    sc.settings.verbosity = 2

    # Ensure correct types for paths and not tuples
    for k in ["path002", "source_he", "sr_spatial"]:
        v = cfg[k]
        if isinstance(v, tuple):
            raise TypeError(f"cfg['{k}'] is a tuple (likely a shift), not a path: {v}")
        if not isinstance(v, (str, Path)):
            raise TypeError(f"cfg['{k}'] must be a str or Path, got {type(v)}: {v}")
            
    # Load + basic filter + destripe
    log_debug(f"Loading sample {sample_key} from {cfg['path002']}")
    adata = b2c.read_visium(
        cfg["path002"],
        source_image_path=cfg["source_he"],
        spaceranger_image_path=cfg["sr_spatial"]
    )
    log_debug(f"Loaded AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")

    adata.var_names_make_unique()
    log_debug("Made var names unique.")

    n_cells_pre = adata.n_obs
    n_genes_pre = adata.n_vars
    # Gene filtering
    
    #sc.pp.filter_cells(adata, min_counts=1) (remove for Nat5)
    log_debug(f"Before gene filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    sc.pp.filter_genes(adata, min_cells=3)
    log_debug(f"After gene filtering (min_cells=3): {adata.n_obs} cells, {adata.n_vars} genes")

    # Add n_counts if not present
    if 'n_counts' not in adata.obs:
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1 if scipy.sparse.issparse(adata.X) else adata.X.sum(axis=1)
        log_debug("Added n_counts to adata.obs.")
    else:
        log_debug("n_counts already present in adata.obs.")

    # Destripe
    log_debug("Calling b2c.destripe...")
    b2c.destripe(adata)
    log_debug("Finished destriping.")

    # manual pixel shift
    dx, dy = cfg["shift"]
    adata.obsm["spatial"][:,0] += dx
    adata.obsm["spatial"][:,1] += dy
    log_debug(f"Applied manual pixel shift: dx={dx}, dy={dy}")

    # QC plot 1 striped counts before/after shift
    if do_qc:
        b = subset(adata, MASK1)
        log_debug(f"After MASK1 subset: {b.n_obs} cells in region r={MASK1['r']}, c={MASK1['c']}")
        if b.n_obs > 0:
            b2c.scaled_he_image(b, mpp=0.5)
            b = b[b.obs['n_counts_adjusted'] > THR_COUNTS]
            log_debug(f"After n_counts_adjusted filter (>{THR_COUNTS}): {b.n_obs} cells")
            sc.pl.spatial(b, color=[None,"n_counts_adjusted"],
                          img_key="0.5_mpp_150_buffer", basis="spatial_cropped_150_buffer",
                          alpha=0.5, cmap="gist_rainbow", save=f"_{sample_key}_qc1.png")
        else:
            log_debug("Warning: No cells in MASK1 region for QC plot 1, skipping.")

    # H&E Tiff for Stardist
    he_tiff = cfg["outdir"]/ "stardist/he.tiff"
    b2c.scaled_he_image(adata, mpp=MICS_PER_PIXEL, save_path=he_tiff)
    log_debug(f"Saved scaled H&E image to {he_tiff}")

    # Stardist on H&E
    he_npz = cfg["outdir"]/ "stardist/he.npz"
    b2c.stardist(he_tiff, he_npz, "2D_versatile_he", prob_thresh=PROB_HE)
    log_debug(f"Ran Stardist on H&E, output: {he_npz}")
    b2c.insert_labels(adata, str(he_npz), basis="spatial",
                      spatial_key="spatial_cropped_150_buffer",
                      mpp=MICS_PER_PIXEL, labels_key="labels_he")
    log_debug("Inserted Stardist labels into adata.obs['labels_he'].")

    # QC plot 2 nuclei labels in mask2
    if do_qc:
        b = subset(adata, MASK2)
        log_debug(f"After MASK2 subset: {b.n_obs} cells in region r={MASK2['r']}, c={MASK2['c']}")
        if b.n_obs > 0:
            b = b[b.obs['labels_he']>0]
            log_debug(f"After labels_he>0 filter: {b.n_obs} cells")
            b.obs['labels_he'] = b.obs['labels_he'].astype(str)
            sc.pl.spatial(b, color=[None,"labels_he"],
                          img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
                          save=f"_{sample_key}_qc2.png")
            crop = b2c.get_crop(b, basis="spatial",
                                spatial_key="spatial_cropped_150_buffer",
                                mpp=MICS_PER_PIXEL)
            rend = b2c.view_stardist_labels(he_tiff, he_npz, crop=crop)
            rend = np.clip(rend, 0, 1)
            plt.imsave(cfg["outdir"]/f"{sample_key}_qc2_overlay.png", rend)
            log_debug(f"Saved QC2 overlay to {cfg['outdir']}/{sample_key}_qc2_overlay.png")
        else:
            log_debug("Warning: No cells in MASK2 region for QC plot 2, skipping.")

    # Expand labels
    b2c.expand_labels(adata, "labels_he", "labels_he_expanded",
                      max_bin_distance=MAX_BIN_DISTANCE)
    log_debug("Expanded labels_he to labels_he_expanded.")

    # Optional GEX Stardist + Salvage
    gex_img = b2c.grid_image(adata, "n_counts_adjusted",
                             mpp=MICS_PER_PIXEL, sigma=5)
    gex_tiff = cfg["outdir"]/ "stardist/gex.tiff"

    if gex_img is not None:
        gex_img_uint8 = (255 * (gex_img - gex_img.min()) / (gex_img.ptp() + 1e-8)).astype(np.uint8)
        cv2.imwrite(str(gex_tiff), gex_img_uint8)
        log_debug(f"Saved GEX tiff to {gex_tiff}")
    else:
        log_debug("Warning: gex_img is None, skipping cv2.imwrite.")

    gex_npz  = cfg["outdir"]/ "stardist/gex.npz"
    b2c.stardist(gex_tiff, gex_npz, "2D_versatile_fluo",
                 prob_thresh=PROB_GEX, nms_thresh=0.1)
    log_debug(f"Ran Stardist on GEX, output: {gex_npz}")
    b2c.insert_labels(adata, str(gex_npz), basis="array",
                      mpp=MICS_PER_PIXEL, labels_key="labels_gex")
    b2c.salvage_secondary_labels(adata, "labels_he_expanded", "labels_gex",
                                 labels_key="labels_joint")
    log_debug("Inserted and salvaged secondary labels.")

    # QC plt 3 showing label source in mask2
    if do_qc:
        b = subset(adata, MASK2)
        log_debug(f"After MASK2 subset for QC3: {b.n_obs} cells")
        if b.n_obs > 0:
            b = b[b.obs["labels_joint"]>0]
            log_debug(f"After labels_joint>0 filter: {b.n_obs} cells")
            b.obs["labels_joint"] = b.obs["labels_joint"].astype(str)
            sc.pl.spatial(b, color=[None,"labels_joint_source","labels_joint"],
                          img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
                          save=f"_{sample_key}_qc3.png")
        else:
            log_debug("Warning: No cells in MASK2 region for QC plot 3, skipping.")

    # Bin to cell
    cdata = b2c.bin_to_cell(adata, "labels_joint",
                            spatial_keys=["spatial","spatial_cropped_150_buffer"])
    log_debug(f"Ran bin_to_cell: {cdata.n_obs} cells, {cdata.n_vars} genes")

    # Round the counts to the nearest integer
    if isinstance(cdata.X, np.ndarray):
        cdata.X = np.round(cdata.X)
        log_debug("Rounded cdata.X (dense array)")
    elif scipy.sparse.issparse(cdata.X):
        cdata.X.data = np.round(cdata.X.data)
        log_debug("Rounded cdata.X (sparse array)")
    else:
        log_debug(f"Warning: Could not round cdata.X: {type(cdata.X)}")

    # QC plot 4 cell-level crop
    if do_qc:
        d = subset(cdata, MASK2)
        log_debug(f"After MASK2 subset for QC4: {d.n_obs} cells")
        if d.n_obs > 0:
            sc.pl.spatial(d, color=["bin_count"],
                          img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
                          save=f"_{sample_key}_qc4.png")
        else:
            log_debug("Warning: No cells in MASK2 region for QC plot 4, skipping.")

    adata.write_h5ad(cfg["outdir"]/f"{sample_key}_2um_raw.h5ad")
    cdata.write_h5ad(cfg["outdir"]/f"{sample_key}_b2c_cells.h5ad")
    log_debug(f"Saved AnnData and cdata to {cfg['outdir']}")
    print(f" {sample_key}: {cdata.n_obs:,} cells – finished")

    # Write debug log to markdown file
    debug_md = cfg["outdir"]/"bin2cell_debug_log.md"
    with open(debug_md, "w") as f:
        f.write("# Bin2Cell Debug Log\n\n")
        for line in debug_logs:
            f.write(f"- {line}\n")
    print(f"Debug log written to {debug_md}")

# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sample", choices=SAMPLES,
                        help="key from the SAMPLES dictionary")
    parser.add_argument("--qc", action="store_true",
                        help="generate the four PNG QC figures")
    args = parser.parse_args()

    run(args.sample, args.qc)
