#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/5_select.out
#SBATCH --error=logs/5_select.out

set -euo pipefail

CONDA_ENV_NAME="${CONDA_ENV_NAME:-ted}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export MPLCONFIGDIR PYTHONUNBUFFERED=1
mkdir -p logs results/plots/select "${MPLCONFIGDIR}"

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

python src/5_select_rank/select_candidates.py --n 25 --max-identity 0.8 --min-iptm 0.5
python src/5_select_rank/summarize_manifest.py
python src/5_select_rank/analyze_binder_selection.py
