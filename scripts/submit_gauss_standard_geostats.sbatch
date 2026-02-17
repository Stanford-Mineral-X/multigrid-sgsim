#!/bin/bash
#SBATCH --job-name=gauss_geostats
#SBATCH --partition=serc
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-32
#SBATCH --output=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/gauss_geostats_%A_%a.out
#SBATCH --error=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/gauss_geostats_%A_%a.err

# ============================================================================
# Standard Geostatistics (Kriging + SGSIM) — GAUSSIAN VARIOGRAM
# ============================================================================
#
# Same structure as submit_standard_geostats.sbatch but uses Gaussian
# variogram configs and saves results to gauss_vgram/ subdirectory.
#
# Array structure (33 tasks):
#   Task 0:       Kriging for dense data
#   Task 1:       Kriging for sparse data
#   Task 2:       Kriging for medium data
#   Tasks 3-12:   SGSIM realizations for dense (100 per task = 1000 total)
#   Tasks 13-22:  SGSIM realizations for sparse (100 per task = 1000 total)
#   Tasks 23-32:  SGSIM realizations for medium (100 per task = 1000 total)
#
# ============================================================================

# --- USER PARAMETERS ---
REALIZATIONS_PER_ENSEMBLE=1000
SGSIM_TASKS_PER_DATASET=10
REALIZATIONS_PER_TASK=$((REALIZATIONS_PER_ENSEMBLE / SGSIM_TASKS_PER_DATASET))
OUTPUT_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results/gauss_vgram/standard_geostats"

# Repo layout (works regardless of where sbatch is submitted from)
REPO_DIR="/oak/stanford/groups/cyaolai/JoshRines/repos/multigrid-sgsim"
SCRIPT_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/configs"
# Gaussian variogram configs
CONFIG_KRIGING_DENSE="${CONFIG_DIR}/config_gauss_standard_kriging_dense.json"
CONFIG_KRIGING_SPARSE="${CONFIG_DIR}/config_gauss_standard_kriging_sparse.json"
CONFIG_KRIGING_MEDIUM="${CONFIG_DIR}/config_gauss_standard_kriging_medium.json"
CONFIG_SGSIM_DENSE="${CONFIG_DIR}/config_gauss_standard_sgsim_dense.json"
CONFIG_SGSIM_SPARSE="${CONFIG_DIR}/config_gauss_standard_sgsim_sparse.json"
CONFIG_SGSIM_MEDIUM="${CONFIG_DIR}/config_gauss_standard_sgsim_medium.json"

# --- ENVIRONMENT SETUP ---
echo "=============================================="
echo "Standard Geostatistics Processing (GAUSSIAN)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo "=============================================="

# Load environment
ml system
ml python/3.12.1
source /oak/stanford/groups/cyaolai/JoshRines/mineralx/venvs/venv_mgsim/bin/activate

# Force unbuffered Python output (for live log viewing)
export PYTHONUNBUFFERED=1

# Create output directories
mkdir -p ${OUTPUT_DIR}/dense
mkdir -p ${OUTPUT_DIR}/sparse
mkdir -p ${OUTPUT_DIR}/medium

# Create log directory for this job
mkdir -p /oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_${SLURM_ARRAY_JOB_ID}

# Change to script directory
cd $SCRIPT_DIR

# --- DETERMINE WHAT TO RUN BASED ON ARRAY TASK ID ---
TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ ${TASK_ID} -eq 0 ]; then
    # Task 0: Kriging for dense
    echo ""
    echo "=== Running Kriging on DENSE data (Gaussian) ==="
    echo "Config: ${CONFIG_KRIGING_DENSE}"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_KRIGING_DENSE} \
        --method kriging \
        --output ${OUTPUT_DIR}/dense/kriging_result.nc

elif [ ${TASK_ID} -eq 1 ]; then
    # Task 1: Kriging for sparse
    echo ""
    echo "=== Running Kriging on SPARSE data (Gaussian) ==="
    echo "Config: ${CONFIG_KRIGING_SPARSE}"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_KRIGING_SPARSE} \
        --method kriging \
        --output ${OUTPUT_DIR}/sparse/kriging_result.nc

elif [ ${TASK_ID} -eq 2 ]; then
    # Task 2: Kriging for medium
    echo ""
    echo "=== Running Kriging on MEDIUM data (Gaussian) ==="
    echo "Config: ${CONFIG_KRIGING_MEDIUM}"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_KRIGING_MEDIUM} \
        --method kriging \
        --output ${OUTPUT_DIR}/medium/kriging_result.nc

elif [ ${TASK_ID} -ge 3 ] && [ ${TASK_ID} -le 12 ]; then
    # Tasks 3-12: SGSIM for dense
    SGSIM_TASK=$((TASK_ID - 3))  # 0-9
    START=$((SGSIM_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${SGSIM_TASK} -eq $((SGSIM_TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running SGSIM on DENSE data (Gaussian) ==="
    echo "Config: ${CONFIG_SGSIM_DENSE}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_SGSIM_DENSE} \
        --method sgsim \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/dense/sgsim_realizations_${START}_${END}.nc

elif [ ${TASK_ID} -ge 13 ] && [ ${TASK_ID} -le 22 ]; then
    # Tasks 13-22: SGSIM for sparse
    SGSIM_TASK=$((TASK_ID - 13))  # 0-9
    START=$((SGSIM_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${SGSIM_TASK} -eq $((SGSIM_TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running SGSIM on SPARSE data (Gaussian) ==="
    echo "Config: ${CONFIG_SGSIM_SPARSE}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_SGSIM_SPARSE} \
        --method sgsim \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/sparse/sgsim_realizations_${START}_${END}.nc

elif [ ${TASK_ID} -ge 23 ] && [ ${TASK_ID} -le 32 ]; then
    # Tasks 23-32: SGSIM for medium
    SGSIM_TASK=$((TASK_ID - 23))  # 0-9
    START=$((SGSIM_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${SGSIM_TASK} -eq $((SGSIM_TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running SGSIM on MEDIUM data (Gaussian) ==="
    echo "Config: ${CONFIG_SGSIM_MEDIUM}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_SGSIM_MEDIUM} \
        --method sgsim \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/medium/sgsim_realizations_${START}_${END}.nc

else
    echo "ERROR: Unknown task ID ${TASK_ID}"
    exit 1
fi

# --- COMPLETION ---
echo ""
echo "=============================================="
echo "Task ${TASK_ID} completed at: $(date)"
echo "=============================================="
