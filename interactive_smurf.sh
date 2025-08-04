#!/bin/bash
#SBATCH -J p5_smurf_jupyter
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=0-10:00:00
#SBATCH --output=interactive_smurfjupyter_%j.log
#SBATCH --error=interactive_smurfjupyter_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jjoseph@nygenome.org

set -euo pipefail

echo "==== Job environment ===="
echo "JobID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Working dir: $PWD"
echo "========================="

module purge
module load miniconda3
module load cuda/11.7.0
module load cudnn/8.5.0.96-CUDA-11.7.0

eval "$(conda shell.bash hook)"
conda activate /gpfs/commons/groups/innovation/jjoseph/enact_results/integrate_env

# Make SMURF importable regardless of env; adjust path if needed
export PYTHONPATH="/gpfs/commons/groups/innovation/jjoseph/smurf_results/SMURF:${PYTHONPATH:-}"

export MPLBACKEND=Agg

PORT=${PORT:-8888}
if ss -tln | grep -q ":$PORT "; then
  PORT=$(shuf -i 8000-9999 -n 1)
fi

echo "Starting Jupyter on port: $PORT"
echo "Connect with (from your laptop):"
echo "  ssh -N -L 8888:$(hostname -s):$PORT $USER@ne1-login01"
echo
echo "Then open: http://localhost:8888/?token=... (see token below)"
echo

which jupyter
jupyter lab --no-browser --ip=0.0.0.0 --port=${PORT} --NotebookApp.allow_remote_access=True
