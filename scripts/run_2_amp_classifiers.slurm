#!/bin/bash
#SBATCH --job-name=classifiers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --output=logs/classifiers.out
#SBATCH --error=logs/classifiers.out
#
# Run AMP Scanner, Macrel, HydrAMP, APEX, and ToxinPred on FASTA files.
#
# Tunables (override at submit time):
#   FASTAS_DIR (default: fastas)          directory of input .fasta files
#   OUTPUT_DIR (default: results/classifiers)
#   ORCHESTRATOR_ENV (default: gmconda_py3923)  conda env for run_classifiers.py
#   CLASSIFIERS (default: ampscanner macrel hydramp_amp hydramp_mic toxinpred)
#
# PREREQUISITE done on a LOGIN NODE (compute nodes usually have no internet):
#   cd src/2_amp_classifiers && ./download_models.sh && ./setup_envs.sh

FASTAS_DIR="${FASTAS_DIR:-fastas}"
OUTPUT_DIR="${OUTPUT_DIR:-results/classifiers}"
# HydrAMP is DROPPED by default. It only accepts peptides of <= 25 aa
# (HYDRAMP_MAX_LEN in merge_classifier_results.py); with the binder range now
# 50-100 aa every sequence is filtered out, the Keras model receives an empty
# batch and dies with "Unexpected result of `predict_function` (Empty
# batch_outputs)", taking the whole stage - and the chain - down with it.
#
# Nothing depends on it: no gate in filters.json uses a hydramp column, and
# ranking.py carries hydramp_* only as passthrough annotations. It is an
# AMP-likeness classifier, which is explicitly not an objective here.
# merge_classifier_results.py handles its absence natively - returns {} for the
# missing CSV and marks every row hydramp_scored=False with
# hydramp_skip_reason="length>25", which satisfies create_manifest's validation.
# Re-add with CLASSIFIERS="... hydramp_amp hydramp_mic ..." if you go back to
# short peptides.
CLASSIFIERS="${CLASSIFIERS:-ampscanner macrel toxinpred}"

set -euo pipefail

REPO_DIR="$(pwd)"
mkdir -p "${REPO_DIR}/${OUTPUT_DIR}" "${REPO_DIR}/logs"

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate gmconda_py3923
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

FASTA_FILES=( "${REPO_DIR}/${FASTAS_DIR}"/*.fasta )

cd "${REPO_DIR}/src/2_amp_classifiers"
for fasta in "${FASTA_FILES[@]}"; do
    echo "=== $(basename "${fasta}") ($(date)) ==="
    ./run_classifiers.sh \
        --fasta "${fasta}" \
        --output-dir "${REPO_DIR}/${OUTPUT_DIR}" \
        --classifiers ${CLASSIFIERS}
done

echo "Merging classifier outputs ..."
python "${REPO_DIR}/src/2_amp_classifiers/merge_classifier_results.py" \
    --results-dir "${REPO_DIR}/${OUTPUT_DIR}" \
    --fastas-dir "${REPO_DIR}/${FASTAS_DIR}" \
    --output "${REPO_DIR}/${OUTPUT_DIR}/merged_classifier_results.csv"

echo "Classifiers finished at $(date). Outputs in ${REPO_DIR}/${OUTPUT_DIR}/"
