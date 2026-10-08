#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/1_design-%j.out
#SBATCH --error=logs/1_design-%j.out

set -euo pipefail

REPO_DIR="$(pwd)"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-boltzgen}"
NUM_DESIGNS="${NUM_DESIGNS:-5000}"
BUDGET="${BUDGET:-200}"
TARGETS="${TARGETS:-KPC3 NDM5}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/scratch/tmp/mpl-${USER:-ampbinder}}"
export MPLCONFIGDIR PYTHONUNBUFFERED=1

mkdir -p "${REPO_DIR}/logs" "${MPLCONFIGDIR}"

echo "========== Stage 1: design =========="
echo "Job ID: ${SLURM_JOB_ID:-N/A} | Node: $(hostname) | Date: $(date)"
echo "Repo: ${REPO_DIR}"
echo "Env: ${CONDA_ENV_NAME} | num_designs=${NUM_DESIGNS} | budget=${BUDGET}"
echo "Targets: ${TARGETS}"
echo "====================================="

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_BASE}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${REPO_DIR}/src/1_design:${PYTHONPATH:-}"

command -v boltzgen >/dev/null 2>&1 || {
  echo "boltzgen CLI not found in conda env ${CONDA_ENV_NAME}" >&2
  exit 1
}

python - <<'PY'
from pathlib import Path
import boltzgen
import boltzgen.cli.boltzgen as cli

paths = [Path(boltzgen.__file__).resolve(), Path(cli.__file__).resolve(), Path(cli.main_script).resolve()]
for path in paths:
    print(f"BoltzGen path: {path}")
    if "/scratch/esevinc22/pepgen/" in str(path):
        raise SystemExit(f"Refusing to run with PepGen-linked BoltzGen path: {path}")
if not all("/AMPBinderDesign/AMPBinderDesign/src/1_design/boltzgen/" in str(path) for path in paths):
    raise SystemExit(f"Unexpected BoltzGen import path(s): {paths}")
PY

# Refuse to design if a target's structure and design sequence disagree. The
# stage-0 scaffolds are NDM-1 (5YPM) and KPC-2 (3DW0); design_domain.fasta
# carries the NDM-5 / KPC-3 mutations but the coordinates did not, so the two
# arms silently targeted different proteins at interface positions.
# Set TARGET_CHECK=warn to downgrade to a warning.
python src/0_targets/check_target_consistency.py --targets ${TARGETS} \
  $([[ "${TARGET_CHECK:-strict}" == "warn" ]] && echo --warn-only) || {
  echo "Target structure/sequence mismatch - run src/0_targets/apply_mutations.py first." >&2
  exit 1
}

# BoltzGen design path.
for target in ${TARGETS}; do
  yaml="boltzgen_inputs/yamls/${target}.yaml"
  out_dir="boltzgen_outputs/${target}"
  if [[ ! -f "${yaml}" ]]; then
    echo "Missing design YAML: ${yaml}" >&2
    exit 1
  fi
  echo "=== BoltzGen: ${target} ($(date)) ==="
  # --additional_filters 'CYS_fraction<0' adds the hard no-designed-cysteine
  # filter. BoltzGen only enables filter_cysteine for the peptide-anything /
  # antibody / nanobody protocols; the default protein-anything protocol we use
  # does NOT, which is why the original run produced panels where 10/25 NDM5
  # designs carried a free cysteine. ('<' means lower-is-better.)
  boltzgen run "${yaml}" \
    --output "${out_dir}" \
    --protocol "${BOLTZGEN_PROTOCOL:-protein-anything}" \
    --num_designs "${NUM_DESIGNS}" \
    --budget "${BUDGET}" \
    --filter_biased true \
    --additional_filters 'CYS_fraction<0' \
    --reuse
done

EXTRACT_OUTPUT_DIRS=(boltzgen_outputs)

# RFdiffusion backbone path. Set RUN_RFDIFFUSION=1 and RFDIFFUSION_REPO to a
# RosettaCommons/RFdiffusion checkout with weights available in models/.
if [[ "${RUN_RFDIFFUSION:-0}" == "1" ]]; then
  : "${RFDIFFUSION_REPO:?Set RFDIFFUSION_REPO=/path/to/RFdiffusion}"
  RFDIFFUSION_ENV="${RFDIFFUSION_ENV:-SE3nv}"
  set +u
  conda activate "${RFDIFFUSION_ENV}"
  set -u
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
  python src/1_design/rfdiffusion.py \
    --rfdiffusion-dir "${RFDIFFUSION_REPO}" \
    --model-dir "${RFDIFFUSION_MODEL_DIR:-${RFDIFFUSION_REPO}/models}" \
    --out-dir rfdiffusion_outputs \
    --num-designs "${RFDIFFUSION_NUM_DESIGNS:-1000}" \
    --binder-len "${RFDIFFUSION_BINDER_LEN:-50-100}" \
    --hotspot-limit "${RFDIFFUSION_HOTSPOT_LIMIT:-0}" \
    --noise "${RFDIFFUSION_NOISE:-0.5}"

  set +u
  conda activate "${CONDA_ENV_NAME}"
  set -u
  export PYTHONPATH="${REPO_DIR}/src/1_design:${PYTHONPATH:-}"
  python src/1_design/inverse_folding.py --design_dir rfdiffusion_outputs/KPC3 --output_dir rfdiffusion_outputs/KPC3_inverse_folding --suffix ".pdb" --skip_existing --multiplicity "${IFOLD_MULTIPLICITY:-10}"
  python src/1_design/inverse_folding.py --design_dir rfdiffusion_outputs/NDM5 --output_dir rfdiffusion_outputs/NDM5_inverse_folding --suffix ".pdb" --skip_existing --multiplicity "${IFOLD_MULTIPLICITY:-10}"
  EXTRACT_OUTPUT_DIRS+=(rfdiffusion_outputs/KPC3_inverse_folding rfdiffusion_outputs/NDM5_inverse_folding)
fi

python src/1_design/extract_binders.py \
  --outputs-dir "${EXTRACT_OUTPUT_DIRS[@]}" \
  --out-dir fastas
