#!/bin/bash
# Stage 1d: put the RFdiffusion designs through BoltzGen's OWN quality control.
#
# Why
# ---
# The BoltzGen arm is filtered by `boltzgen run`'s filter step: it refolds every
# design in complex, design-chain-only, and in isolation, applies hard
# thresholds, then ranks. It rejects ~90-96% of its own output. The RFdiffusion
# arm had no quality control at all, so it entered the pool unselected and
# outnumbered the QC-passed BoltzGen arm ~20-50x.
#
# From the `folding` step onward BoltzGen's pipeline operates on a single
# directory of inverse-folded designs (<id>.cif + <id>.npz), and
# rfdiffusion_outputs/<T>_inverse_folding/ already has exactly that layout with
# identical .npz keys. So we can run steps folding/design_folding/analysis/
# filtering over it and get an all_designs_metrics.csv in the same format,
# which extract_binders.py then gates on (QC_METRICS in that file).
#
# The design CIF/NPZ files are HARDLINKED into a separate tree so QC artifacts
# (refold_cif/, fold_out_npz/, ...) do not pollute the design output dir.
#
# Env:
#   TARGETS            default "KPC3 NDM5"
#   RFQC_BUDGET        diversity-selected set size (default 500)
#   RFQC_LIMIT         only QC the first N designs per target (0 = all).
#                      Use this for a pilot before committing GPU-days.
#   BOLTZGEN_PROTOCOL  default protein-anything (needed for design_folding)
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/1d_rfdiffusion_qc-%j.out
#SBATCH --error=logs/1d_rfdiffusion_qc-%j.out

set -euo pipefail

REPO_DIR="$(pwd)"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltzgen}"
TARGETS="${TARGETS:-KPC3 NDM5}"
RFQC_BUDGET="${RFQC_BUDGET:-500}"
RFQC_LIMIT="${RFQC_LIMIT:-0}"
BOLTZGEN_PROTOCOL="${BOLTZGEN_PROTOCOL:-protein-anything}"
NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/scratch/tmp/numba-cache-${USER:-ampbinder}}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export NUMBA_CACHE_DIR MPLCONFIGDIR PYTHONUNBUFFERED=1
mkdir -p logs "${NUMBA_CACHE_DIR}" "${MPLCONFIGDIR}"

echo "========== Stage 1d: BoltzGen QC on RFdiffusion designs =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Targets: ${TARGETS} | budget=${RFQC_BUDGET} | limit=${RFQC_LIMIT:-all} | protocol=${BOLTZGEN_PROTOCOL}"
echo "=================================================================="

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${REPO_DIR}/src/1_design:${PYTHONPATH:-}"

command -v boltzgen >/dev/null 2>&1 || { echo "boltzgen CLI not found" >&2; exit 1; }

for target in ${TARGETS}; do
  src_dir="rfdiffusion_outputs/${target}_inverse_folding"
  qc_root="rfdiffusion_qc/${target}"
  dest_dir="${qc_root}/intermediate_designs_inverse_folded"
  yaml="boltzgen_inputs/yamls/${target}.yaml"

  [[ -d "${src_dir}" ]] || { echo "Missing ${src_dir}; run stage 1b first" >&2; exit 1; }
  [[ -f "${yaml}"   ]] || { echo "Missing ${yaml}" >&2; exit 1; }

  mkdir -p "${dest_dir}"
  echo "=== ${target}: staging designs into ${dest_dir} ==="
  n=0
  shopt -s nullglob
  _cifs=( "${src_dir}"/*.cif )
  shopt -u nullglob
  for cif in "${_cifs[@]}"; do
    stem="$(basename "${cif}" .cif)"
    [[ -f "${src_dir}/${stem}.npz" ]] || continue
    if [[ "${RFQC_LIMIT}" != "0" && "${n}" -ge "${RFQC_LIMIT}" ]]; then break; fi
    ln -f "${cif}"                  "${dest_dir}/${stem}.cif" 2>/dev/null || cp -p "${cif}" "${dest_dir}/${stem}.cif"
    ln -f "${src_dir}/${stem}.npz"  "${dest_dir}/${stem}.npz" 2>/dev/null || cp -p "${src_dir}/${stem}.npz" "${dest_dir}/${stem}.npz"
    n=$((n+1))
  done
  echo "    staged ${n} designs (hardlinked, no extra disk)"

  echo "=== ${target}: boltzgen QC steps (folding, design_folding, analysis, filtering) ==="
  boltzgen run "${yaml}" \
    --output "${qc_root}" \
    --steps folding design_folding analysis filtering \
    --protocol "${BOLTZGEN_PROTOCOL}" \
    --budget "${RFQC_BUDGET}" \
    --filter_biased true \
    --additional_filters 'CYS_fraction<0' \
    --reuse

  metrics="${qc_root}/final_ranked_designs/all_designs_metrics.csv"
  if [[ -f "${metrics}" ]]; then
    python - "${metrics}" "${target}" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
ok = sum(1 for r in rows
         if str(r.get("pass_filters", "")).strip().lower() in ("true", "1", "1.0"))
print(f"[{sys.argv[2]}] RFdiffusion QC: {ok}/{len(rows)} designs passed "
      f"({100*ok/max(len(rows),1):.1f}%)")
PY
  else
    echo "[${target}] WARNING: ${metrics} was not written - QC did not complete" >&2
  fi
done

echo "=== done ($(date)). extract_binders.py will now gate rfdiffusion_ifold on these verdicts. ==="
