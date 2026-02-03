#!/bin/bash
#SBATCH --job-name=std_geostats
#SBATCH --partition=serc
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-9
#SBATCH --output=logs/std_geostats_%A_%a.out
#SBATCH --error=logs/std_geostats_%A_%a.err

# =============================================================================
# Run Standard Geostatistics (Kriging + SGSIM) for Comparison
#
# This script runs as a SLURM array job:
# - Task 0: Kriging for all datasets
# - Tasks 1-9: SGSIM realizations (100 per task = 1000 total per dataset)
#
# Usage:
#   sbatch submit_standard_geostats.sbatch
# =============================================================================

# --- Configuration ---
DATASETS=("dense" "sparse")
REALIZATIONS_PER_TASK=100
CONFIGS_DIR="../configs"
RESULTS_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results/standard_geostats"

# --- Environment ---
module load python/3.9
source ~/envs/mgsim/bin/activate  # Adjust to your conda/venv

# --- Create directories ---
mkdir -p ${RESULTS_DIR}
mkdir -p logs

# --- Determine what to run based on array task ID ---
TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ ${TASK_ID} -eq 0 ]; then
    # Task 0: Run kriging for all datasets
    echo "=== Running Kriging for all datasets ==="

    for DATASET in "${DATASETS[@]}"; do
        CONFIG="${CONFIGS_DIR}/config_iso_subregions_${DATASET}.json"
        OUTPUT_DIR="${RESULTS_DIR}/${DATASET}"
        mkdir -p ${OUTPUT_DIR}

        echo "Running Kriging on ${DATASET}..."
        python run_standard_geostats.py \
            --config ${CONFIG} \
            --method kriging \
            --output ${OUTPUT_DIR}/kriging_result.nc

        echo "Kriging ${DATASET} complete."
    done

else
    # Tasks 1-9: Run SGSIM realizations
    SGSIM_TASK=$((TASK_ID - 1))  # 0-8
    START=$((SGSIM_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    echo "=== Running SGSIM realizations ${START}-${END} ==="

    for DATASET in "${DATASETS[@]}"; do
        CONFIG="${CONFIGS_DIR}/config_iso_subregions_${DATASET}.json"
        OUTPUT_DIR="${RESULTS_DIR}/${DATASET}"
        mkdir -p ${OUTPUT_DIR}

        echo "Running SGSIM ${START}-${END} on ${DATASET}..."
        python run_standard_geostats.py \
            --config ${CONFIG} \
            --method sgsim \
            --start ${START} \
            --end ${END} \
            --output ${OUTPUT_DIR}/sgsim_realizations_${START}_${END}.nc

        echo "SGSIM ${DATASET} ${START}-${END} complete."
    done
fi

echo "Task ${TASK_ID} finished at $(date)"
