#!/bin/bash
# Stage 0: make the target structures match the design sequence, then
# regenerate the BoltzGen design specs.
#
# 1. apply_mutations.py  - the stage-0 scaffolds are NDM-1 (5YPM) and KPC-2
#    (3DW0). design_domain.fasta carries the NDM-5 / KPC-3 mutations but the
#    coordinates never did, and prep_report.md referenced an apply_mutations.py
#    that did not exist. The affected residues sit ON the designed interface
#    (NDM5 V88 3.6 A, M154 3.7 A, KPC3 H274 is itself an epitope residue), so
#    the RFdiffusion arm was designing against the wrong protein.
#    NOTE: new sidechains are truncated to backbone+CB - no coordinates are
#    invented. Rebuild with a packer if you need full sidechains.
# 2. check_target_consistency.py - hard gate; the design scripts run it too.
# 3. make_config.py - regenerates boltzgen_inputs/yamls/<T>.yaml so the new
#    BINDER_LEN actually reaches BoltzGen. This writer previously emitted
#    boltzgen_inputs/<T>_design.yaml, which nothing read.
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --output=logs/0_targets-%j.out
#SBATCH --error=logs/0_targets-%j.out

set -euo pipefail

CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltzgen}"
TARGETS="${TARGETS:-KPC3 NDM5}"
export PYTHONUNBUFFERED=1
mkdir -p logs

echo "========== Stage 0: target structures + design specs =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Targets: ${TARGETS}"
echo "==============================================================="

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

echo "--- 1. applying documented mutations to the structures ---"
for target in ${TARGETS}; do
  python src/0_targets/apply_mutations.py \
    --target "${target}" \
    --accept-truncated \
    --out "targets/${target}/structures/${target}_target_mutated.pdb"
done

echo "--- 2. structure / sequence consistency ---"
python src/0_targets/check_target_consistency.py --targets ${TARGETS}

echo "--- 3. regenerating BoltzGen design specs (picks up BINDER_LEN) ---"
python src/1_design/make_config.py
for target in ${TARGETS}; do
  echo "    ${target} binder spec: $(grep -A2 'id: P' boltzgen_inputs/yamls/${target}.yaml | grep 'sequence:' | tr -s ' ')"
done

echo "=== stage 0 done ($(date)) ==="
