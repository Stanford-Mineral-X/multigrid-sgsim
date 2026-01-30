#!/bin/bash
#SBATCH --job-name=mgsim
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-9
#SBATCH --output=logs/mgsim_%A_%a.out
#SBATCH --error=logs/mgsim_%A_%a.err

# ============================================================================
# MGSIM Batch Processing on Sherlock
# ============================================================================
#
# This script runs MGSIM realizations using SLURM array jobs.
# Each array task handles a subset of realizations.
#
# Usage:
#   1. Edit the parameters below
#   2. Create logs directory: mkdir -p logs
#   3. Submit: sbatch submit_mgsim.sbatch
#
# To run 1000 realizations with 10 array tasks:
#   - Set TOTAL_REALIZATIONS=1000
#   - Set #SBATCH --array=0-9
#   - Each task will run 100 realizations
#
# ============================================================================

# --- USER PARAMETERS (edit these) ---
TOTAL_REALIZATIONS=1000
CONFIG_FILE="mgsim_config.pkl"
OUTPUT_DIR="results"

# --- COMPUTED PARAMETERS ---
NUM_TASKS=10  # Should match --array range (0-9 = 10 tasks)
REALIZATIONS_PER_TASK=$((TOTAL_REALIZATIONS / NUM_TASKS))

START_IDX=$((SLURM_ARRAY_TASK_ID * REALIZATIONS_PER_TASK))
END_IDX=$((START_IDX + REALIZATIONS_PER_TASK))

# Handle last task (might have extra realizations due to rounding)
if [ $SLURM_ARRAY_TASK_ID -eq $((NUM_TASKS - 1)) ]; then
    END_IDX=$TOTAL_REALIZATIONS
fi

OUTPUT_FILE="${OUTPUT_DIR}/mgsim_realizations_${START_IDX}_${END_IDX}.nc"

# --- ENVIRONMENT SETUP ---
echo "=============================================="
echo "MGSIM Batch Processing"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""
echo "Processing realizations $START_IDX to $((END_IDX - 1))"
echo "Output file: $OUTPUT_FILE"
echo "=============================================="

# Load conda/mamba environment
# Adjust this to match your Sherlock setup
module load anaconda3
source activate mgsim  # or: conda activate mgsim

# Create output directory
mkdir -p $OUTPUT_DIR
mkdir -p logs

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
