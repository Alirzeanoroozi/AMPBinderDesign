#!/bin/bash
#SBATCH --job-name=run_1b_rfdiffusion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=logs/1b_rfdiffusion-%j.out
#SBATCH --error=logs/1b_rfdiffusion-%j.out

set -euo pipefail

REPO_DIR="$(pwd)"
RFDIFFUSION_REPO="${RFDIFFUSION_REPO:-/scratch/esevinc22/RFdiffusion}"
RFDIFFUSION_MODEL_DIR="${RFDIFFUSION_MODEL_DIR:-${RFDIFFUSION_REPO}/models}"
RFDIFFUSION_ENV="${RFDIFFUSION_ENV:-SE3nv}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltzgen}"
TARGETS="${TARGETS:-KPC3 NDM5}"
RFDIFFUSION_NUM_DESIGNS="${RFDIFFUSION_NUM_DESIGNS:-1000}"
RFDIFFUSION_BINDER_LEN="${RFDIFFUSION_BINDER_LEN:-50-100}"
RFDIFFUSION_HOTSPOT_LIMIT="${RFDIFFUSION_HOTSPOT_LIMIT:-0}"
RFDIFFUSION_NOISE="${RFDIFFUSION_NOISE:-0.5}"
MOL_DIR="${MOL_DIR:-/home/esevinc22/.cache/huggingface/hub/datasets--boltzgen--inference-data/snapshots/c3d36fd276e9caf098c75d4113c6d5eb320b1a4c/mols.zip}"
IFOLD_CHECKPOINT="${IFOLD_CHECKPOINT:-/home/esevinc22/.cache/huggingface/hub/models--boltzgen--boltzgen-1/snapshots/c1be29e1f82ffcc72264f64b993c43fb4e0d17f0/boltzgen1_ifold.ckpt}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export MPLCONFIGDIR PYTHONUNBUFFERED=1

mkdir -p "${REPO_DIR}/logs" "${MPLCONFIGDIR}"

echo "========== Stage 1b: RFDiffusion =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Repo: ${REPO_DIR}"
echo "RFdiffusion: ${RFDIFFUSION_REPO}"
echo "Targets: ${TARGETS} | num_designs=${RFDIFFUSION_NUM_DESIGNS}"
echo "==========================================="


CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"

set +u
conda activate "${RFDIFFUSION_ENV}"
set -u

# Refuse to design if a target's structure and design sequence disagree. The
# stage-0 scaffolds are NDM-1 (5YPM) and KPC-2 (3DW0); design_domain.fasta
# carries the NDM-5 / KPC-3 mutations but the coordinates did not, so this arm
# silently designed against the wrong residue at interface positions. Must run
# AFTER conda activation - `python` does not exist before it.
# Set TARGET_CHECK=warn to downgrade to a warning.
rc=0
python src/0_targets/check_target_consistency.py --targets ${TARGETS} \
  $([[ "${TARGET_CHECK:-strict}" == "warn" ]] && echo --warn-only) || rc=$?
if [[ "${rc}" -ne 0 ]]; then
  echo "Target consistency check failed (exit ${rc}; 1=mismatch, 2=cannot check)." >&2
  echo "Run: python src/0_targets/apply_mutations.py --target <T> --accept-truncated \\" >&2
  echo "       --out targets/<T>/structures/<T>_target_mutated.pdb" >&2
  exit 1
fi
python - <<'PY'
import importlib.util
import sys
import torch

required = ["rfdiffusion", "dgl", "hydra", "omegaconf", "se3_transformer"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    raise SystemExit(f"Missing RFdiffusion dependencies: {', '.join(missing)}")
if torch.version.cuda is None:
    raise SystemExit("RFdiffusion env has CPU-only Torch; install a CUDA-enabled Torch build.")
if not torch.cuda.is_available():
    raise SystemExit("CUDA Torch is installed, but no CUDA device is visible to this job.")
print(
    "RFdiffusion env OK:",
    f"python={sys.executable}",
    f"torch={torch.__version__}",
    f"cuda={torch.version.cuda}",
    f"device={torch.cuda.get_device_name(0)}",
)
PY
GEN_TARGETS=""
for target in ${TARGETS}; do
  # nullglob array, not `ls | wc -l`: on a FRESH run the directory does not
  # exist yet, ls exits 2, pipefail propagates it and set -e kills the job
  # silently (2>/dev/null hides the reason).
  shopt -s nullglob
  _raw=( "rfdiffusion_outputs/${target}/${target}"_*.pdb )
  shopt -u nullglob
  existing=${#_raw[@]}
  if [[ "${existing}" -lt "${RFDIFFUSION_NUM_DESIGNS}" ]]; then
    GEN_TARGETS="${GEN_TARGETS} ${target}"
  else
    echo "Skipping RFdiffusion generation for ${target}: ${existing}/${RFDIFFUSION_NUM_DESIGNS} raw designs already on disk"
  fi
done

if [[ -n "${GEN_TARGETS}" ]]; then
  python src/1_design/rfdiffusion.py \
    --targets ${GEN_TARGETS} \
    --rfdiffusion-dir "${RFDIFFUSION_REPO}" \
    --model-dir "${RFDIFFUSION_MODEL_DIR}" \
    --out-dir rfdiffusion_outputs \
    --num-designs "${RFDIFFUSION_NUM_DESIGNS}" \
    --binder-len "${RFDIFFUSION_BINDER_LEN}" \
    --hotspot-limit "${RFDIFFUSION_HOTSPOT_LIMIT}" \
    --noise "${RFDIFFUSION_NOISE}"
fi

set +u
conda activate "${CONDA_ENV_NAME}"
set -u
export PYTHONPATH="${REPO_DIR}/src/1_design:${PYTHONPATH:-}"
IFOLD_MULTIPLICITY="${IFOLD_MULTIPLICITY:-10}"
for target in ${TARGETS}; do
  design_dir="rfdiffusion_outputs/${target}"
  output_dir="rfdiffusion_outputs/${target}_inverse_folding"
  pending_dir="rfdiffusion_outputs/${target}_ifold_pending"

  # inverse_folding.py's own --skip_existing compares against "<stem>.cif",
  # but this pipeline writes one file per sequence variant, "<stem>_N.cif"
  # for N in [0, multiplicity), so it never matches. A design is only truly
  # done once all IFOLD_MULTIPLICITY variant files exist; pre-filter to a
  # symlink dir of designs missing any variant so a rerun resumes correctly.
  rm -rf "${pending_dir}"
  mkdir -p "${pending_dir}"
  pending_count=0
  # nullglob throughout: an unguarded glob iterates the LITERAL pattern when it
  # matches nothing, which would yield a bogus stem like "KPC3_*" and then
  # symlink a non-existent file into the pending dir.
  shopt -s nullglob
  _designs=( "${design_dir}/${target}"_*.pdb )
  for pdb in "${_designs[@]}"; do
    stem="$(basename "${pdb}" .pdb)"
    _variants=( "${output_dir}/${stem}_"*.cif )
    have=${#_variants[@]}
    if [[ "${have}" -lt "${IFOLD_MULTIPLICITY}" ]]; then
      _parts=( "${design_dir}/${stem}".* )
      for f in "${_parts[@]}"; do
        ln -sf "$(readlink -f "${f}")" "${pending_dir}/"
      done
      pending_count=$((pending_count + 1))
    fi
  done
  shopt -u nullglob
  if [[ "${#_designs[@]}" -eq 0 ]]; then
    echo "${target}: no RFdiffusion backbones found in ${design_dir}" >&2
    exit 1
  fi
  echo "${target}: ${pending_count}/${RFDIFFUSION_NUM_DESIGNS} designs still need inverse folding (need ${IFOLD_MULTIPLICITY} variants each)"

  if [[ "${pending_count}" -gt 0 ]]; then
    # --multiplicity MUST be passed: inverse_folding.py defaults to 10, so
    # IFOLD_MULTIPLICITY previously only affected the "is this design done?"
    # bookkeeping above and not generation. Setting it to 5 silently still
    # produced 10 variants per backbone; setting it above 10 would have made the
    # pending loop never converge.
    python src/1_design/inverse_folding.py \
      --design_dir "${pending_dir}" \
      --output_dir "${output_dir}" \
      --moldir "${MOL_DIR}" \
      --checkpoint "${IFOLD_CHECKPOINT}" \
      --multiplicity "${IFOLD_MULTIPLICITY}" \
      --suffix ".pdb"
  else
    echo "${target}: inverse folding already complete, skipping"
  fi
  rm -rf "${pending_dir}"
done
