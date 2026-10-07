#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/4_boltz2.out
#SBATCH --error=logs/4_boltz2.out

set -euo pipefail

CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltz}"
NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/scratch/tmp/numba-cache-${USER:-ampbinder}}"
BOLTZ_CACHE="${BOLTZ_CACHE:-/scratch/esevinc22/.boltz}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export NUMBA_CACHE_DIR BOLTZ_CACHE MPLCONFIGDIR PYTHONUNBUFFERED=1
mkdir -p logs results "${NUMBA_CACHE_DIR}" "${BOLTZ_CACHE}" "${MPLCONFIGDIR}"

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

BOLTZ_DEVICES="${BOLTZ_DEVICES:-1}"
# Where to read filtered_<TARGET>.fasta from. Point at a subset directory
# (e.g. results/stage4_sample) to predict a stratified sample instead of the
# whole stage-3 survivor pool.
FILTERED_DIR="${FILTERED_DIR:-results}"

python src/4_structure_prediction/prepare_configs.py --filtered-dir "${FILTERED_DIR}"
boltz predict boltz_inputs --out_dir boltz_results --write_full_pae --devices "${BOLTZ_DEVICES}"

PRED=boltz_results/boltz_results_boltz_inputs/predictions

python src/4_structure_prediction/parse_boltz2.py --boltz_out boltz_results --target NDM5
python src/4_structure_prediction/parse_boltz2.py --boltz_out boltz_results --target KPC3

python src/4_structure_prediction/active_site_overlap.py --target NDM5 --complexes "${PRED}"
python src/4_structure_prediction/active_site_overlap.py --target KPC3 --complexes "${PRED}"

python src/4_structure_prediction/run_ipsae.py --boltz-out boltz_results --workers 8

python src/4_structure_prediction/build_manifest.py
