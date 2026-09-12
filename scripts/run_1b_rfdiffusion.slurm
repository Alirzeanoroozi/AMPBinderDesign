#!/bin/bash
#SBATCH --job-name=run_1b_rfdiffusion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=200G
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
RFDIFFUSION_BINDER_LEN="${RFDIFFUSION_BINDER_LEN:-12-45}"
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

conda activate "${RFDIFFUSION_ENV}"
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
  --targets ${TARGETS} \
  --rfdiffusion-dir "${RFDIFFUSION_REPO}" \
  --model-dir "${RFDIFFUSION_MODEL_DIR}" \
  --out-dir rfdiffusion_outputs \
  --num-designs "${RFDIFFUSION_NUM_DESIGNS}" \
  --binder-len "${RFDIFFUSION_BINDER_LEN}" \
  --hotspot-limit "${RFDIFFUSION_HOTSPOT_LIMIT}" \
  --noise "${RFDIFFUSION_NOISE}"

conda activate "${CONDA_ENV_NAME}"
export PYTHONPATH="${REPO_DIR}/src/1_design:${PYTHONPATH:-}"
for target in ${TARGETS}; do
  python src/1_design/inverse_folding.py \
    --design_dir "rfdiffusion_outputs/${target}" \
    --output_dir "rfdiffusion_outputs/${target}_inverse_folding" \
    --moldir "${MOL_DIR}" \
    --checkpoint "${IFOLD_CHECKPOINT}" \
    --suffix ".pdb" \
    --skip_existing
done
