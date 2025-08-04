#!/usr/bin/env python3

#file for generating configs with specific parametesr that you can run on batch script, after this run a batch script to run them all 

import yaml
import os

# Base paths
BASE_DATA_PATH = "/gpfs/commons/groups/innovation/jjoseph/data"
BASE_RESULTS_PATH = "/gpfs/commons/groups/innovation/jjoseph/enact_results"

#list samples you want to run 
# Sample configurations
SAMPLES = [
    {"name": "P2_CRC", "analysis_name": "colon_P2CRC_gpu", "type": "Cancer"},
    {"name": "P5_CRC", "analysis_name": "colon_P5CRC_gpu", "type": "Cancer"},
]

# Bin-to-cell assignment methods
BIN_TO_CELL_METHODS = ["naive", "weighted_by_area", "weighted_by_gene", "weighted_by_cluster"]

NBINS_LIST = [2, 4, 6, 8, 10]
DEFAULT_NBINS = 2  # Default nbins for the second set

# Load base config as template
with open("config/configs.yaml", "r") as f:
    base_config = yaml.safe_load(f)

# weighted_by_area with different nbins for both P2 and P5
print("=== Generating Set 1: weighted_by_area with different nbins ===")
for sample in SAMPLES:
    for nbins in NBINS_LIST:
        config = base_config.copy()
        config["analysis_name"] = f"{sample['analysis_name']}_weighted_by_area_nbins{nbins}"
        config["cache_dir"] = os.path.join(BASE_RESULTS_PATH, sample["name"], "Cache")
        sample_prefix = f"Visium_HD_Human_Colon_{sample['type']}_{sample['name'].split('_')[0]}" 
        config["paths"]["wsi_path"] = os.path.join(BASE_DATA_PATH, sample["name"], "inputs", f"{sample_prefix}_tissue_image.btf")
        config["paths"]["visiumhd_h5_path"] = os.path.join(BASE_DATA_PATH, sample["name"], "binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
        config["paths"]["tissue_positions_path"] = os.path.join(BASE_DATA_PATH, sample["name"], "binned_outputs/square_002um/spatial/tissue_positions.parquet")
        config["params"]["bin_to_cell_method"] = "weighted_by_area"
        config["params"]["expand_by_nbins"] = nbins
        output_path = f"config/configs_{sample['name']}_weighted_by_area_nbins{nbins}.yaml"
        with open(output_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
        print(f"Generated config for {sample['name']} with weighted_by_area and nbins {nbins}: {output_path}")

# All methods with default nbins for both P2 and P5
print("\n=== Generating Set 2: All methods with default nbins ===")
for sample in SAMPLES:
    for method in BIN_TO_CELL_METHODS:
        config = base_config.copy()
        config["analysis_name"] = f"{sample['analysis_name']}_{method}_default"
        config["cache_dir"] = os.path.join(BASE_RESULTS_PATH, sample["name"], "Cache")
        sample_prefix = f"Visium_HD_Human_Colon_{sample['type']}_{sample['name'].split('_')[0]}" 
        config["paths"]["wsi_path"] = os.path.join(BASE_DATA_PATH, sample["name"], "inputs", f"{sample_prefix}_tissue_image.btf")
        config["paths"]["visiumhd_h5_path"] = os.path.join(BASE_DATA_PATH, sample["name"], "binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
        config["paths"]["tissue_positions_path"] = os.path.join(BASE_DATA_PATH, sample["name"], "binned_outputs/square_002um/spatial/tissue_positions.parquet")
        config["params"]["bin_to_cell_method"] = method
        config["params"]["expand_by_nbins"] = DEFAULT_NBINS
        output_path = f"config/configs_{sample['name']}_{method}_default.yaml"
        with open(output_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
        print(f"Generated config for {sample['name']} with {method} and default nbins: {output_path}")
        