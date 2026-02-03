#!/bin/bash
#SBATCH --job-name=std_geostats
#SBATCH --partition=serc
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-21
#SBATCH --output=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/std_geostats_%A_%a.out
#SBATCH --error=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/std_geostats_%A_%a.err

# ============================================================================
# Standard Geostatistics (Kriging + SGSIM) for Comparison with MGSIM
# ============================================================================
#
# This script runs Kriging and SGSIM using the SAME config files as MGSIM
# to ensure a fair comparison (same data, same grid, same variograms).
#
# Array structure (22 tasks):
#   Task 0:       Kriging for dense data
#   Task 1:       Kriging for sparse data
#   Tasks 2-11:   SGSIM realizations for dense (100 per task = 1000 total)
#   Tasks 12-21:  SGSIM realizations for sparse (100 per task = 1000 total)
#
# We use the iso_subregions configs since Kriging/SGSIM use a single global
# variogram anyway (not cluster-specific).
#
# ============================================================================

# --- USER PARAMETERS ---
REALIZATIONS_PER_ENSEMBLE=1000
SGSIM_TASKS_PER_DATASET=10
REALIZATIONS_PER_TASK=$((REALIZATIONS_PER_ENSEMBLE / SGSIM_TASKS_PER_DATASET))
OUTPUT_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results/standard_geostats"

# Config files (same as MGSIM ensembles for data consistency)
CONFIG_DIR="${SLURM_SUBMIT_DIR}/../configs"
CONFIG_DENSE="${CONFIG_DIR}/config_iso_subregions_dense.json"
CONFIG_SPARSE="${CONFIG_DIR}/config_iso_subregions_sparse.json"

# --- ENVIRONMENT SETUP ---
echo "=============================================="
echo "Standard Geostatistics Processing"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo "=============================================="

# Load environment (same as submit_all_ensembles.sbatch)
ml system
ml python/3.12.1
source /oak/stanford/groups/cyaolai/JoshRines/mineralx/venvs/venv_mgsim/bin/activate

# Force unbuffered Python output (for live log viewing)
export PYTHONUNBUFFERED=1

# Create output directories
mkdir -p ${OUTPUT_DIR}/dense
mkdir -p ${OUTPUT_DIR}/sparse

# Create log directory for this job
mkdir -p /oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_${SLURM_ARRAY_JOB_ID}

# Change to script directory
cd $SLURM_SUBMIT_DIR

# --- DETERMINE WHAT TO RUN BASED ON ARRAY TASK ID ---
TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ ${TASK_ID} -eq 0 ]; then
    # Task 0: Kriging for dense
    echo ""
    echo "=== Running Kriging on DENSE data ==="
    echo "Config: ${CONFIG_DENSE}"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_DENSE} \
        --method kriging \
        --output ${OUTPUT_DIR}/dense/kriging_result.nc

elif [ ${TASK_ID} -eq 1 ]; then
    # Task 1: Kriging for sparse
    echo ""
    echo "=== Running Kriging on SPARSE data ==="
    echo "Config: ${CONFIG_SPARSE}"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_SPARSE} \
        --method kriging \
        --output ${OUTPUT_DIR}/sparse/kriging_result.nc

elif [ ${TASK_ID} -ge 2 ] && [ ${TASK_ID} -le 11 ]; then
    # Tasks 2-11: SGSIM for dense
    SGSIM_TASK=$((TASK_ID - 2))  # 0-9
    START=$((SGSIM_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${SGSIM_TASK} -eq $((SGSIM_TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running SGSIM on DENSE data ==="
    echo "Config: ${CONFIG_DENSE}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_DENSE} \
        --method sgsim \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/dense/sgsim_realizations_${START}_${END}.nc

elif [ ${TASK_ID} -ge 12 ] && [ ${TASK_ID} -le 21 ]; then
    # Tasks 12-21: SGSIM for sparse
    SGSIM_TASK=$((TASK_ID - 12))  # 0-9
    START=$((SGSIM_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${SGSIM_TASK} -eq $((SGSIM_TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running SGSIM on SPARSE data ==="
    echo "Config: ${CONFIG_SPARSE}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_SPARSE} \
        --method sgsim \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/sparse/sgsim_realizations_${START}_${END}.nc

else
    echo "ERROR: Unknown task ID ${TASK_ID}"
    exit 1
fi

# --- COMPLETION ---
echo ""
echo "=============================================="
echo "Task ${TASK_ID} completed at: $(date)"
echo "=============================================="
