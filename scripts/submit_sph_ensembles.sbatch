#!/bin/bash
#SBATCH --job-name=sph_mgsim
#SBATCH --partition=serc
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-49
#SBATCH --output=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/sph_mgsim_%A_%a.out
#SBATCH --error=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/sph_mgsim_%A_%a.err

# ============================================================================
# MGSIM 5-Ensemble Batch Processing — SPHERICAL VARIOGRAM
# ============================================================================
#
# Same structure as submit_all_ensembles.sbatch but uses spherical variogram
# configs (config_sph_*.json) and saves results to sph_vgram/ subdirectory.
#
# Array structure (50 tasks):
#   Array tasks 0-9:   sph_iso_subregions_dense   (100 realizations each)
#   Array tasks 10-19: sph_iso_global_dense       (100 realizations each)
#   Array tasks 20-29: sph_aniso_subregions_dense (100 realizations each)
#   Array tasks 30-39: sph_aniso_global_dense     (100 realizations each)
#   Array tasks 40-49: sph_iso_subregions_sparse  (100 realizations each)
#
# Total: 50 array tasks x 100 = 5000 realizations (1000 per ensemble)
#
# ============================================================================

# --- USER PARAMETERS ---
REALIZATIONS_PER_ENSEMBLE=1000
TASKS_PER_ENSEMBLE=10
OUTPUT_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results/sph_vgram"

# Repo layout (works regardless of where sbatch is submitted from)
REPO_DIR="/oak/stanford/groups/cyaolai/JoshRines/repos/multigrid-sgsim"
SCRIPT_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/configs"

# --- COMPUTED PARAMETERS ---
REALIZATIONS_PER_TASK=$((REALIZATIONS_PER_ENSEMBLE / TASKS_PER_ENSEMBLE))

# Determine which ensemble and which task within that ensemble
ENSEMBLE_IDX=$((SLURM_ARRAY_TASK_ID / TASKS_PER_ENSEMBLE))
TASK_WITHIN_ENSEMBLE=$((SLURM_ARRAY_TASK_ID % TASKS_PER_ENSEMBLE))

# Map ensemble index to name (spherical configs use sph_ prefix)
case $ENSEMBLE_IDX in
    0) ENSEMBLE_NAME="iso_subregions_dense" ;;
    1) ENSEMBLE_NAME="iso_global_dense" ;;
    2) ENSEMBLE_NAME="aniso_subregions_dense" ;;
    3) ENSEMBLE_NAME="aniso_global_dense" ;;
    4) ENSEMBLE_NAME="iso_subregions_sparse" ;;
esac

CONFIG_FILE="${CONFIG_DIR}/config_sph_${ENSEMBLE_NAME}.json"

# Calculate realization range for this task
START_IDX=$((TASK_WITHIN_ENSEMBLE * REALIZATIONS_PER_TASK))
END_IDX=$((START_IDX + REALIZATIONS_PER_TASK))

# Handle last task in ensemble
if [ $TASK_WITHIN_ENSEMBLE -eq $((TASKS_PER_ENSEMBLE - 1)) ]; then
    END_IDX=$REALIZATIONS_PER_ENSEMBLE
fi

# Output file
ENSEMBLE_OUTPUT_DIR="${OUTPUT_DIR}/${ENSEMBLE_NAME}"
OUTPUT_FILE="${ENSEMBLE_OUTPUT_DIR}/realizations_${START_IDX}_${END_IDX}.nc"

# --- ENVIRONMENT SETUP ---
echo "=============================================="
echo "MGSIM Ensemble Processing (SPHERICAL VARIOGRAM)"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""
echo "Ensemble: sph_$ENSEMBLE_NAME"
echo "Config: $CONFIG_FILE"
echo "Realizations: $START_IDX to $((END_IDX - 1))"
echo "Output: $OUTPUT_FILE"
echo "=============================================="

# Load environment
ml system
ml python/3.12.1
source /oak/stanford/groups/cyaolai/JoshRines/mineralx/venvs/venv_mgsim/bin/activate

# Force unbuffered Python output (for live log viewing)
export PYTHONUNBUFFERED=1

# Create output directories
mkdir -p $ENSEMBLE_OUTPUT_DIR

# Create log directory for this job
mkdir -p /oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_${SLURM_ARRAY_JOB_ID}

# Change to script directory
cd $SCRIPT_DIR

# --- RUN MGSIM ---
python run_mgsim_batch.py \
    --config $CONFIG_FILE \
    --start $START_IDX \
    --end $END_IDX \
    --output $OUTPUT_FILE

# --- COMPLETION ---
echo ""
echo "=============================================="
echo "Completed at: $(date)"
echo "=============================================="
