#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/boltzgen.out
#SBATCH --error=logs/boltzgen.out

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate bg_macrel
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

# BoltzGen design path.
# python src/1_design/design.py --yaml_path boltzgen_inputs/yamls/KPC3.yaml --output_dir boltzgen_outputs/KPC3/
# python src/1_design/design.py --yaml_path boltzgen_inputs/yamls/NDM5.yaml --output_dir boltzgen_outputs/NDM5/

# python src/1_design/inverse_folding.py --design_dir boltzgen_outputs/KPC3 --output_dir boltzgen_outputs/KPC3_inverse_folding/
# python src/1_design/inverse_folding.py --design_dir boltzgen_outputs/NDM5 --output_dir boltzgen_outputs/NDM5_inverse_folding/

EXTRACT_OUTPUT_DIRS=(boltzgen_outputs)

# RFdiffusion backbone path. Set RUN_RFDIFFUSION=1 and RFDIFFUSION_REPO to a
# RosettaCommons/RFdiffusion checkout with weights available in models/.
if [[ "${RUN_RFDIFFUSION:-0}" == "1" ]]; then
  : "${RFDIFFUSION_REPO:?Set RFDIFFUSION_REPO=/path/to/RFdiffusion}"
  RFDIFFUSION_ENV="${RFDIFFUSION_ENV:-SE3nv}"
  source ${CONDA_BASE}/bin/activate "${RFDIFFUSION_ENV}"
  python src/1_design/rfdiffusion.py \
    --rfdiffusion-dir "${RFDIFFUSION_REPO}" \
    --model-dir "${RFDIFFUSION_MODEL_DIR:-${RFDIFFUSION_REPO}/models}" \
    --out-dir rfdiffusion_outputs \
    --num-designs "${RFDIFFUSION_NUM_DESIGNS:-1000}" \
    --binder-len "${RFDIFFUSION_BINDER_LEN:-12-45}" \
    --hotspot-limit "${RFDIFFUSION_HOTSPOT_LIMIT:-0}" \
    --noise "${RFDIFFUSION_NOISE:-0.5}"

  source ${CONDA_BASE}/bin/activate bg_macrel
  python src/1_design/inverse_folding.py --design_dir rfdiffusion_outputs/KPC3 --output_dir rfdiffusion_outputs/KPC3_inverse_folding --skip_existing
  python src/1_design/inverse_folding.py --design_dir rfdiffusion_outputs/NDM5 --output_dir rfdiffusion_outputs/NDM5_inverse_folding --skip_existing
  EXTRACT_OUTPUT_DIRS+=(rfdiffusion_outputs/KPC3_inverse_folding rfdiffusion_outputs/NDM5_inverse_folding)
fi

python src/1_design/extract_binders.py \
  --outputs-dir "${EXTRACT_OUTPUT_DIRS[@]}" \
  --out-dir fastas
