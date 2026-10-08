#!/bin/bash
# Rebuild the structure manifest (now joining interface_metrics_<T>.csv) and
# re-select with the multimetric peptide-aware ranking. CPU only - no GPU, no
# new Boltz-2 inference.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/5b_rebuild_select-%j.out
#SBATCH --error=logs/5b_rebuild_select-%j.out

set -euo pipefail
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export MPLCONFIGDIR PYTHONUNBUFFERED=1
mkdir -p logs results/plots/select "${MPLCONFIGDIR}"

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME:-ted}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

python src/4_structure_prediction/build_manifest.py
python src/5_select_rank/select_candidates.py --n 25 --max-identity 0.8 --min-iptm 0.5
python src/5_select_rank/summarize_manifest.py
python src/5_select_rank/analyze_binder_selection.py
