#!/bin/bash
#SBATCH --job-name=csb_mgsim
#SBATCH --partition=serc
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-49
#SBATCH --output=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/csb_mgsim_%A_%a.out
#SBATCH --error=/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/logs/job_%A/csb_mgsim_%A_%a.err

# ============================================================================
# MGSIM Batch Processing — Cape Smith Belt (CSB) Real Dataset
# ============================================================================
#
# Runs 1000 realizations of MGSIM with:
#   - Spherical isotropic variograms (10 clusters / subregions)
#   - RBF trend surface (pre-computed in config)
#   - mg_resols = [500, 300, 250, 200, 150, 100, 50]
#   - num_points = 20, radius = 500
#
# Array structure (50 tasks):
#   Each task runs 20 realizations sequentially.
#   Task 0:  realizations 0-19
#   Task 1:  realizations 20-39
#   ...
#   Task 49: realizations 980-999
#
# Estimated time: ~25 min/realization × 20 = ~8.3 hours per task
# Total: 50 array tasks × 20 = 1000 realizations
#
# ============================================================================

# --- USER PARAMETERS ---
REALIZATIONS_TOTAL=1000
TASKS_TOTAL=50
CONFIG_NAME="csb_sph_iso"
OUTPUT_DIR="/oak/stanford/groups/cyaolai/JoshRines/sherlock/sherlock_multigrid/results/csb_sph_iso"

# Repo layout
REPO_DIR="/oak/stanford/groups/cyaolai/JoshRines/repos/multigrid-sgsim"
SCRIPT_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/configs"

# --- COMPUTED PARAMETERS ---
REALIZATIONS_PER_TASK=$((REALIZATIONS_TOTAL / TASKS_TOTAL))

CONFIG_FILE="${CONFIG_DIR}/config_${CONFIG_NAME}.json"

# Calculate realization range for this task
START_IDX=$((SLURM_ARRAY_TASK_ID * REALIZATIONS_PER_TASK))
END_IDX=$((START_IDX + REALIZATIONS_PER_TASK))

# Handle last task (pick up any remainder)
if [ $SLURM_ARRAY_TASK_ID -eq $((TASKS_TOTAL - 1)) ]; then
    END_IDX=$REALIZATIONS_TOTAL
fi

# Output file
OUTPUT_FILE="${OUTPUT_DIR}/realizations_${START_IDX}_${END_IDX}.nc"

# --- ENVIRONMENT SETUP ---
echo "=============================================="
echo "MGSIM — Cape Smith Belt (CSB) Real Dataset"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""
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

# Create output directory
mkdir -p $OUTPUT_DIR

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
