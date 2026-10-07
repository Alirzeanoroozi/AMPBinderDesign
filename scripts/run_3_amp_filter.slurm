#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/3_amp_filter.out
#SBATCH --error=logs/3_amp_filter.out

set -euo pipefail

CONDA_ENV_NAME="${CONDA_ENV_NAME:-ted}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export MPLCONFIGDIR PYTHONUNBUFFERED=1
mkdir -p logs results/plots "${MPLCONFIGDIR}"

# Initialise environment and modules
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

python src/3_amp_filter/developability.py --fasta fastas/NDM5.fasta --out results/developability_NDM5.csv
python src/3_amp_filter/developability.py --fasta fastas/KPC3.fasta --out results/developability_KPC3.csv

python src/3_amp_filter/create_manifest.py \
  --classifiers results/classifiers/merged_classifier_results.csv

python src/3_amp_filter/plot_distribution.py

python src/3_amp_filter/apply_filters.py --manifest results/manifest_NDM5.csv --out results/filtered_NDM5.csv
python src/3_amp_filter/apply_filters.py --manifest results/manifest_KPC3.csv --out results/filtered_KPC3.csv
