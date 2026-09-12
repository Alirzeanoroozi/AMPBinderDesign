#!/bin/bash
#SBATCH --job-name=run_1c_extract_fastas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/1c_extract_fastas-%j.out
#SBATCH --error=logs/1c_extract_fastas-%j.out

set -euo pipefail

REPO_DIR="$(pwd)"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltzgen}"
TARGETS="${TARGETS:-KPC3 NDM5}"
EXTRACT_OUTPUT_DIRS=(boltzgen_outputs)

for target in ${TARGETS}; do
  EXTRACT_OUTPUT_DIRS+=("rfdiffusion_outputs/${target}_inverse_folding")
done

mkdir -p "${REPO_DIR}/logs" "${REPO_DIR}/fastas"

echo "========== Stage 1c: extract FASTAs =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Repo: ${REPO_DIR}"
echo "Input roots: ${EXTRACT_OUTPUT_DIRS[*]}"
echo "=============================================="

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

python src/1_design/extract_binders.py \
  --outputs-dir "${EXTRACT_OUTPUT_DIRS[@]}" \
  --out-dir fastas \
  --targets ${TARGETS}
