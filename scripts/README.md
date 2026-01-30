# MGSIM HPC Batch Processing Scripts

Scripts for running MGSIM realizations on Sherlock (Stanford HPC) or any SLURM-based cluster.

## Overview

| Script | Purpose |
|--------|---------|
| `prepare_mgsim_config.py` | Prepares input data and saves config pickle (single ensemble) |
| `prepare_ensemble_configs.py` | Prepares 4 config files for ensemble comparison |
| `run_mgsim_batch.py` | Runs batch realizations (called by SLURM) |
| `submit_mgsim.sbatch` | SLURM batch script for single ensemble |
| `submit_all_ensembles.sbatch` | SLURM batch script for all 4 ensembles |
| `combine_results.py` | Combines partial NetCDF files into one |
| `combine_and_compare_ensembles.py` | Combines all ensembles and generates comparison |
| `compute_error_metrics.py` | Computes error metrics vs ground truth |

## 4-Ensemble Comparison (Recommended)

For the paper, run 4 ensembles to compare:
1. **Isotropic + Subregions**: Isotropic variograms, cluster-specific
2. **Isotropic + Global**: Isotropic variograms, single global
3. **Anisotropic + Subregions**: Directional variograms, cluster-specific
4. **Anisotropic + Global**: Directional variograms, single global

### Quick Start (4 Ensembles)

```bash
# 1. Prepare all 4 configs
python prepare_ensemble_configs.py

# 2. On Sherlock: Submit all 4 ensembles (4000 total realizations)
mkdir -p logs results
sbatch submit_all_ensembles.sbatch

# 3. After completion: Combine and compare
python combine_and_compare_ensembles.py \
    --results-dir results \
    --ground-truth /path/to/gt_xyvc.csv \
    --output-dir analysis
```

This generates:
- Combined NetCDF for each ensemble
- Error maps comparing all 4
- Box plots of metric distributions
- LaTeX table for the paper (`comparison_table.tex`)

## Quick Start

### 1. Prepare Configuration (on local machine or login node)

Edit `prepare_mgsim_config.py` to match your dataset, then run:

```bash
python prepare_mgsim_config.py
```

This creates `mgsim_config.pkl` containing:
- Input DataFrame (`df_xyvtcs`)
- Variogram parameters (`df_gamma`)
- Multigrid resolutions
- Grid coordinates
- NST transformer (if using)

### 2. Transfer to Sherlock

```bash
scp -r scripts/ sherlock:/path/to/your/project/
scp mgsim_config.pkl sherlock:/path/to/your/project/scripts/
```

### 3. Submit Jobs

On Sherlock:

```bash
cd /path/to/your/project/scripts
mkdir -p logs results

# Edit submit_mgsim.sbatch to set:
#   - TOTAL_REALIZATIONS (e.g., 1000)
#   - Array size (e.g., --array=0-9 for 10 tasks)
#   - Partition, time limit, memory

sbatch submit_mgsim.sbatch
```

### 4. Monitor Jobs

```bash
squeue -u $USER               # Check job status
tail -f logs/mgsim_*.out      # Watch output
```

### 5. Combine Results

After all jobs complete:

```bash
python combine_results.py --input-dir results --output mgsim_combined.nc
```

### 6. Compute Error Metrics (for synthetic data)

```bash
python compute_error_metrics.py \
    --results mgsim_combined.nc \
    --ground-truth /path/to/gt_xyvc.csv \
    --output error_analysis.nc \
    --figures figures
```

## Output Format

### Combined NetCDF (`mgsim_combined.nc`)

Dimensions:
- `realization`: realization index (0 to N-1)
- `y`: y coordinates
- `x`: x coordinates

Variables:
- `simulated(realization, y, x)`: All realizations
- `mean(y, x)`: Mean across realizations
- `variance(y, x)`: Variance across realizations
- `std(y, x)`: Standard deviation
- `median(y, x)`: Median (50th percentile)
- `p05, p25, p75, p95(y, x)`: Percentiles
- `iqr(y, x)`: Interquartile range

### Error Analysis NetCDF (`error_analysis.nc`)

Variables:
- `ground_truth(y, x)`: Ground truth values
- `mean_prediction(y, x)`: MGSIM mean prediction
- `error(y, x)`: Prediction - Ground truth
- `abs_error(y, x)`: |Error|
- `prediction_variance(y, x)`: Variance from realizations
- `realization_rmse(realization)`: RMSE per realization
- `realization_r2(realization)`: R² per realization
- etc.

## Customization

### Changing Number of Realizations

In `submit_mgsim.sbatch`:

```bash
TOTAL_REALIZATIONS=1000   # Total number of realizations
#SBATCH --array=0-9       # 10 array tasks → 100 realizations each
```

For 1000 realizations with 20 array tasks (50 each):
```bash
TOTAL_REALIZATIONS=1000
#SBATCH --array=0-19
NUM_TASKS=20
```

### Memory/Time Requirements

Adjust based on grid size:
- Small grid (100x100): 4GB, 1 hour per 100 realizations
- Medium grid (500x500): 8GB, 2-4 hours per 100 realizations
- Large grid (1000x1000): 16GB+, may need longer time

### Using NST (Normal Score Transform)

In `prepare_mgsim_config.py`:
```python
use_nst = True
clip_nst = True
clip_percentile = 99.0
```

## Troubleshooting

### Out of Memory
- Increase `#SBATCH --mem`
- Reduce realizations per task (increase array size)

### Jobs Taking Too Long
- Increase `#SBATCH --time`
- Use fewer realizations per task

### Module Not Found
Check conda environment setup in `submit_mgsim.sbatch`:
```bash
module load anaconda3
source activate mgsim
```

## Example Workflow

```bash
# 1. Local: Prepare config
python prepare_mgsim_config.py

# 2. Transfer to Sherlock
rsync -av scripts/ sherlock:~/mgsim/scripts/

# 3. On Sherlock: Submit 1000 realizations (10 array tasks)
ssh sherlock
cd ~/mgsim/scripts
mkdir -p logs results
sbatch submit_mgsim.sbatch

# 4. Wait for completion...
squeue -u $USER

# 5. Combine results
python combine_results.py --input-dir results --output mgsim_1000.nc

# 6. Compute error metrics
python compute_error_metrics.py --results mgsim_1000.nc --ground-truth ../data/gt_xyvc.csv

# 7. Transfer results back
exit
rsync -av sherlock:~/mgsim/scripts/mgsim_1000.nc ./
rsync -av sherlock:~/mgsim/scripts/figures/ ./figures/
```

## 4-Ensemble Workflow (Full Example)

```bash
# 1. Local: Generate 4 config files
python prepare_ensemble_configs.py
# Creates: config_iso_subregions.pkl
#          config_iso_global.pkl
#          config_aniso_subregions.pkl
#          config_aniso_global.pkl

# 2. Transfer to Sherlock
rsync -av scripts/ sherlock:~/mgsim/scripts/

# 3. On Sherlock: Submit all ensembles
ssh sherlock
cd ~/mgsim/scripts
mkdir -p logs results/{iso_subregions,iso_global,aniso_subregions,aniso_global}
sbatch submit_all_ensembles.sbatch
# Submits 40 array tasks: 10 per ensemble × 4 ensembles
# Each task runs 100 realizations → 1000 per ensemble → 4000 total

# 4. Monitor progress
squeue -u $USER
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS

# 5. After completion: Combine and analyze
python combine_and_compare_ensembles.py \
    --results-dir results \
    --ground-truth ../data/gt_xyvc.csv \
    --output-dir analysis

# 6. Transfer results
exit
rsync -av sherlock:~/mgsim/scripts/analysis/ ./analysis/
```

### Output Structure

```
analysis/
├── figures/
│   ├── ensemble_error_maps.png         # Error maps for all 4
│   ├── ensemble_metric_distributions.png  # Box plots
│   └── ensemble_metric_comparison.png  # Bar chart
├── iso_subregions_combined.nc
├── iso_global_combined.nc
├── aniso_subregions_combined.nc
├── aniso_global_combined.nc
├── comparison_table.tex                # LaTeX table for paper
└── ensemble_results.pkl                # Python pickle for further analysis
```

### Expected Results

The comparison should show:
- **Subregions vs Global**: Cluster-specific variograms should perform better in heterogeneous fields
- **Anisotropic vs Isotropic**: Directional variograms should help if the field has preferential orientation

Best case: Anisotropic + Subregions should have lowest RMSE and highest R² if the field exhibits both clustering and anisotropy.
