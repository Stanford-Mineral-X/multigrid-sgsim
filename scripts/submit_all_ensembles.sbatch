#!/bin/bash
#SBATCH --job-name=mgsim_ens
#SBATCH --partition=serc
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-49
#SBATCH --output=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/mgsim_ens_%A_%a.out
#SBATCH --error=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/mgsim_ens_%A_%a.err

# ============================================================================
# MGSIM 5-Ensemble Batch Processing
# ============================================================================
#
# This script runs 5 ensembles x N realizations each:
#   - iso_subregions_dense: Isotropic + cluster-specific variograms + dense FL
#   - iso_global_dense: Isotropic + single global variogram + dense FL
#   - aniso_subregions_dense: Anisotropic + cluster-specific variograms + dense FL
#   - aniso_global_dense: Anisotropic + single global variogram + dense FL
#   - iso_subregions_sparse: Isotropic + cluster-specific + sparse FL (comparison)
#
# Array structure:
#   Array tasks 0-9:   iso_subregions_dense   (100 realizations each)
#   Array tasks 10-19: iso_global_dense       (100 realizations each)
#   Array tasks 20-29: aniso_subregions_dense (100 realizations each)
#   Array tasks 30-39: aniso_global_dense     (100 realizations each)
#   Array tasks 40-49: iso_subregions_sparse  (100 realizations each)
#
# Total: 50 array tasks x 100 = 5000 realizations (1000 per ensemble)
#
# ============================================================================

# --- USER PARAMETERS ---
REALIZATIONS_PER_ENSEMBLE=1000
TASKS_PER_ENSEMBLE=10
OUTPUT_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results"

# --- COMPUTED PARAMETERS ---
REALIZATIONS_PER_TASK=$((REALIZATIONS_PER_ENSEMBLE / TASKS_PER_ENSEMBLE))

# Determine which ensemble and which task within that ensemble
ENSEMBLE_IDX=$((SLURM_ARRAY_TASK_ID / TASKS_PER_ENSEMBLE))
TASK_WITHIN_ENSEMBLE=$((SLURM_ARRAY_TASK_ID % TASKS_PER_ENSEMBLE))

# Map ensemble index to name
case $ENSEMBLE_IDX in
    0) ENSEMBLE_NAME="iso_subregions_dense" ;;
    1) ENSEMBLE_NAME="iso_global_dense" ;;
    2) ENSEMBLE_NAME="aniso_subregions_dense" ;;
    3) ENSEMBLE_NAME="aniso_global_dense" ;;
    4) ENSEMBLE_NAME="iso_subregions_sparse" ;;
esac

CONFIG_FILE="config_${ENSEMBLE_NAME}.pkl"

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
echo "MGSIM Ensemble Processing"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""
echo "Ensemble: $ENSEMBLE_NAME"
echo "Config: $CONFIG_FILE"
echo "Realizations: $START_IDX to $((END_IDX - 1))"
echo "Output: $OUTPUT_FILE"
echo "=============================================="

# Load environment
module load anaconda3
source activate mgsim

# Force unbuffered Python output (for live log viewing)
export PYTHONUNBUFFERED=1

# Create output directories
mkdir -p $ENSEMBLE_OUTPUT_DIR

# Change to script directory
cd $SLURM_SUBMIT_DIR

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
