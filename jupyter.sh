#!/bin/bash 

module load miniconda3
echo "loaded conda"

module load cuda/11.7.0
echo "loaded cuda"

module load cudnn/8.5.0.96-CUDA-11.7.0
echo "loaded cudnn" 

export PYTHONNOUSERSITE=1  # prevents Python from using ~/.local
unset PYTHONPATH           # ensure no local packages override conda env


source activate /gpfs/commons/groups/innovation/jjoseph/enact_results/integrate_env
echo "env activated"


jupyter lab --no-browser --port=8888
echo "set up jupyter"

