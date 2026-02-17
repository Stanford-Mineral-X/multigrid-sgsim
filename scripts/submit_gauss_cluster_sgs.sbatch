#!/bin/bash
#SBATCH --job-name=gauss_csgs
#SBATCH --partition=serc
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-29
#SBATCH --output=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/gauss_csgs_%A_%a.out
#SBATCH --error=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/gauss_csgs_%A_%a.err

# ============================================================================
# Cluster-Specific SGS (No Multigrid) — GAUSSIAN VARIOGRAM
# ============================================================================
#
# Runs SGS with per-cluster variograms in a single pass — no multigrid
# resolution refinement and no trend/residual decomposition.
#
# Workflow (same as standard SGSIM, but per-cluster):
#   1. NST on raw observation values
#   2. Fit a separate variogram per cluster on the NST'd data
#   3. Run cluster_sgs (single pass over full grid)
#   4. Inverse NST to get back to original units
#
# This isolates the effect of cluster-specific variograms vs. the multigrid
# cascade. Variograms are fit at runtime from NST-transformed data per
# cluster (NOT reused from MGSIM, which fits on residuals).
#
# Array structure (30 tasks):
#   Tasks 0-9:   Cluster SGS for dense  (100 realizations each = 1000 total)
#   Tasks 10-19: Cluster SGS for sparse (100 realizations each = 1000 total)
#   Tasks 20-29: Cluster SGS for medium (100 realizations each = 1000 total)
#
# ============================================================================

# --- USER PARAMETERS ---
REALIZATIONS_PER_ENSEMBLE=1000
TASKS_PER_DATASET=10
REALIZATIONS_PER_TASK=$((REALIZATIONS_PER_ENSEMBLE / TASKS_PER_DATASET))
OUTPUT_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results/gauss_vgram/cluster_sgs"

# Repo layout (works regardless of where sbatch is submitted from)
REPO_DIR="/oak/stanford/groups/cyaolai/JoshRines/repos/multigrid-sgsim"
SCRIPT_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/configs"

# Reuse Gaussian SGSIM configs (same data files, variogram model/maxlag/n_lags)
# The cluster_sgs method reads the cluster column from the data CSV and fits
# per-cluster variograms at runtime on NST-transformed values.
CONFIG_DENSE="${CONFIG_DIR}/config_gauss_standard_sgsim_dense.json"
CONFIG_SPARSE="${CONFIG_DIR}/config_gauss_standard_sgsim_sparse.json"
CONFIG_MEDIUM="${CONFIG_DIR}/config_gauss_standard_sgsim_medium.json"

# --- ENVIRONMENT SETUP ---
echo "=============================================="
echo "Cluster-Specific SGS Processing (GAUSSIAN)"
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

if [ ${TASK_ID} -ge 0 ] && [ ${TASK_ID} -le 9 ]; then
    # Tasks 0-9: Cluster SGS for dense
    CSGS_TASK=$((TASK_ID - 0))  # 0-9
    START=$((CSGS_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${CSGS_TASK} -eq $((TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running Cluster SGS on DENSE data (Gaussian) ==="
    echo "Config: ${CONFIG_DENSE}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_DENSE} \
        --method cluster_sgs \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/dense/cluster_sgs_realizations_${START}_${END}.nc

elif [ ${TASK_ID} -ge 10 ] && [ ${TASK_ID} -le 19 ]; then
    # Tasks 10-19: Cluster SGS for sparse
    CSGS_TASK=$((TASK_ID - 10))  # 0-9
    START=$((CSGS_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${CSGS_TASK} -eq $((TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running Cluster SGS on SPARSE data (Gaussian) ==="
    echo "Config: ${CONFIG_SPARSE}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_SPARSE} \
        --method cluster_sgs \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/sparse/cluster_sgs_realizations_${START}_${END}.nc

elif [ ${TASK_ID} -ge 20 ] && [ ${TASK_ID} -le 29 ]; then
    # Tasks 20-29: Cluster SGS for medium
    CSGS_TASK=$((TASK_ID - 20))  # 0-9
    START=$((CSGS_TASK * REALIZATIONS_PER_TASK))
    END=$((START + REALIZATIONS_PER_TASK))

    # Handle last task
    if [ ${CSGS_TASK} -eq $((TASKS_PER_DATASET - 1)) ]; then
        END=${REALIZATIONS_PER_ENSEMBLE}
    fi

    echo ""
    echo "=== Running Cluster SGS on MEDIUM data (Gaussian) ==="
    echo "Config: ${CONFIG_MEDIUM}"
    echo "Realizations: ${START} to $((END - 1))"
    echo ""

    python run_standard_geostats.py \
        --config ${CONFIG_MEDIUM} \
        --method cluster_sgs \
        --start ${START} \
        --end ${END} \
        --output ${OUTPUT_DIR}/medium/cluster_sgs_realizations_${START}_${END}.nc

else
    echo "ERROR: Unknown task ID ${TASK_ID}"
    exit 1
fi

# --- COMPLETION ---
echo ""
echo "=============================================="
echo "Task ${TASK_ID} completed at: $(date)"
echo "=============================================="
