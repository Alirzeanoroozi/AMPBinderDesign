#!/bin/bash
# Stage 4b: peptide-aware interface metrics. CPU ONLY - post-processes the
# CIF/PAE files Boltz-2 already wrote. No GPU, no new inference.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/4b_interface_metrics-%j.out
#SBATCH --error=logs/4b_interface_metrics-%j.out

set -euo pipefail

CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltz}"
export PYTHONUNBUFFERED=1
mkdir -p logs results

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

python src/4_structure_prediction/interface_metrics.py \
  --boltz-out boltz_results \
  --out-dir results \
  --workers "${SLURM_CPUS_PER_TASK:-16}"
