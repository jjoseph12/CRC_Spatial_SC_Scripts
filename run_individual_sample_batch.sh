#!/bin/bash
#SBATCH --job-name=ENACT_${1}_${2}_${3}
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs_07_14/stdout_${1}_${2}_${3}_%j.log
#SBATCH --error=logs_07_14/stderr_${1}_${2}_${3}_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jjoseph@nygenome.org


'''
This script is for running the enact in batch, it works with the generate_config.py although you dont need it and the submit_all_jobs.sh
This file just runs the enact pipeline 
You can clasify here the sample, method and nbins

this is custom script written when doing the qc analysis testing instead of runniing each time

# to run this: ./script.sh <sample> <method> <nbins>

# example: ./script.sh P5 CRC_weighted_by_cluster default

but if you wanted to mass run it, once you have configs set, run the submit_all_jobs.sh file. 

If you dont wnat to use this script which is totally expected, because this was for custom case, 

just pay attention to the exports, and loading modules. Enact ran into some issue with my cluster, and in order to resolve it 
I had to load and export these to get it working. 
Will need to adjust to yoru environment. 
'''


if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <sample> <method> <nbins>"
  exit 1
fi

# Create logs directory
mkdir -p logs_07_14

# Load GPU modules
module purge
module load miniconda3
module load CUDA/12.4.0
module load cuDNN/9.7.0.66-CUDA-12.3.0

export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CUDA_HOME/nvvm/libdevice"
export PYTHONPATH=$HOME/.keras:$PYTHONPATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH


# Most up to date src enact pipeline.py
export PYTHONPATH=/gpfs/commons/groups/innovation/jjoseph/enact-pipeline/src:$PYTHONPATH

# Activate environment
source activate enact_py_env

# Navigate to working directory
cd /gpfs/commons/groups/innovation/jjoseph/enact-pipeline

# XLA workaround
ln -sf $CUDA_HOME/nvvm/libdevice/libdevice.10.bc ./libdevice.10.bc

# Logging
echo ">>> Running ENACT"
echo "Sample: $1"
echo "Method: $2"
echo "Nbins:  $3"

# Config selection logic
if [[ "$2" == "weighted_by_area" && "$3" != "default" ]]; then
  CONFIG_FILE=config/configs_${1}_weighted_by_area_${3}.yaml
elif [[ "$3" == "default" ]]; then
  CONFIG_FILE=config/configs_${1}_${2}_default.yaml
else
  CONFIG_FILE=config/configs_${1}_${2}.yaml
fi

if [[ ! -f $CONFIG_FILE ]]; then
  echo "Config file $CONFIG_FILE does not exist!"
  exit 1
fi

# Run ENACT pipeline
make run_enact CONFIG=$CONFIG_FILE
