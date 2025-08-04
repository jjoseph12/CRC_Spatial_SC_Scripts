#!/bin/bash

#SBATCH --job-name=ENACT_multi_sample           # Job name
#SBATCH --partition=cpu                     # Use CPU partition (no GPU needed)
#SBATCH --mem=64G                          # Request 64GB memory per job
#SBATCH --cpus-per-task=8                  # Request 8 CPU cores per job
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jjoseph@nygenome.org       # Where to send mail  
#SBATCH --time=1-00:00:00                      # Time limit days-hours:minutes:seconds
#SBATCH --output=stdout_multi_%j.log           # output and error log



#this is for generating the config files with CPU, it 
# Load modules 
module load miniconda3

# Load conda environment
source activate enact_py_env

# Change to the pipeline directory
cd /gpfs/commons/groups/innovation/jjoseph/enact-pipeline

# Generate sample-specific configs
python generate_configs.py

# Process each sample
SAMPLES=("P2CRC" "P5CRC")
METHODS=("weighted_by_area" "naive" "weighted_by_cluster" "weighted_by_gene")

for sample in "${SAMPLES[@]}"; do
    for method in "${METHODS[@]}"; do
        echo "Processing sample: $sample with method: $method"
        cp config/configs_${sample}_${method}.yaml config/configs.yaml
        make run_enact
        echo "Completed: $sample with $method"
        echo "----------------------------------------"
    done
done

echo "All runs finished"




