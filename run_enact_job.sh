#!/bin/bash

#SBATCH --job-name=P5_CRC_gpu                 # Job name
#SBATCH --partition=gpu                    # Partition Name
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jjoseph@nygenome.org      # Where to send mail  
#SBATCH --mem=100G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T] 
#Sbatch --time=2-00:00:00                   # Time 2 days
#SBATCH --output=stdout_%j.log  


# This is code uses generic config file config.yml and runs it with with enact platform. 

#load modules 
module load miniconda3

#load the conda environments
source activate enact_py_env

cd /gpfs/commons/groups/innovation/jjoseph/enact-pipeline
make run_enact
