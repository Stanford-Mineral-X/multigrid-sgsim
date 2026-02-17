#!/usr/bin/env python
"""
Prepare MGSIM Configuration for Cape Smith Belt (CSB) Real Dataset — V2

Changes from v1:
  - MAXLAG: 400 -> 2000 (properly capture sill)
  - SMOOTHING: 1000 -> 10000 (smoother trend, larger residuals)
  - SPACING_FL: 200 -> 1000 (sparser trend sampling)

Usage:
    python prepare_csb_config_v2.py
"""

import json
import numpy as np
import pandas as pd
import sys
from pathlib import Path
from scipy.spatial import cKDTree

# Add src directory to path
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

from variograms import cluster_variogram, build_variogram_dataframe
import trendmaking


# =============================================================================
# USER PARAMETERS
# =============================================================================

# Data path
DATA_DIR = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data')
CSB_PATH = DATA_DIR / 'csb_xyvmmc.csv'

# Grid parameters (from notebook: rows=248, cols=426)
ROWS = 248
COLS = 426

# Multigrid resolutions (coarse to fine)
MG_RESOLS = [500, 300, 250, 200, 150, 100, 50]

# MGSIM parameters
NUM_POINTS = 20
RADIUS = 500

# Trend parameters (RBF interpolation) — V2: smoother trend
SMOOTHING = 10000.0
SPACING_FL = 1000

# Variogram fitting parameters — V2: maxlag increased from 400 to 2000
MAXLAG = 2000
N_LAGS = 30
VARIOGRAM_MODEL = 'spherical'

# Output directory for config files
OUTPUT_DIR = script_dir.parent / 'configs'
ENSEMBLE_NAME = 'csb_sph_iso_v2'

# Random seed for trend computation (reproducibility)
RANDOM_SEED = 42


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("CSB MGSIM CONFIGURATION GENERATOR — V2")
    print("  maxlag=2000, smoothing=10000, spacing_fl=1000")
    print("=" * 70)

    # Set random seed for reproducibility of trend
    np.random.seed(RANDOM_SEED)

    # =====================================================================
    # 1. LOAD DATA
    # =====================================================================
    print("\n1. Loading CSB data...")
    df_xyvmmc = pd.read_csv(CSB_PATH)

    # Extract flight lines (fl_mask == 1)
    df_fl_only = df_xyvmmc[df_xyvmmc["fl_mask"] == 1].drop(
        columns=["fl_mask", "grid_mask"], errors="ignore"
    )
    fl_xyvc = df_fl_only.values

    # Extract grid points (drop mask columns)
    df_grid_xyvc = df_xyvmmc.drop(columns=["fl_mask", "grid_mask"])
    grid_xyvc = df_grid_xyvc.values
    grid_xyc = np.column_stack([grid_xyvc[:, :2], grid_xyvc[:, 3]])

    n_regions = len(np.unique(fl_xyvc[:, 3]))
    print(f"  Flight line points: {len(fl_xyvc)}")
    print(f"  Grid points: {len(grid_xyvc)}")
    print(f"  Clusters: {n_regions}")

    # =====================================================================
    # 2. COMPUTE TREND
    # =====================================================================
    print("\n2. Computing RBF trend...")
    fl_xyvct, grid_xyct = trendmaking.make_trend(
        fl_xyvc, grid_xyc, SMOOTHING, SPACING_FL
    )

    # Compute residuals on flight lines
    residuals = fl_xyvct[:, 2] - fl_xyvct[:, 4]
    fl_xyvctr = np.column_stack([fl_xyvct, residuals])

    df_fl_xyvctr = pd.DataFrame(
        fl_xyvctr,
        columns=["x", "y", "value", "cluster", "trend", "residual"]
    )

    df_grid_xyct = pd.DataFrame(
        grid_xyct,
        columns=["x", "y", "cluster", "trend"]
    )

    print(f"  Trend computed. Mean residual: {df_fl_xyvctr['residual'].mean():.3f}")

    # =====================================================================
    # 3. BUILD df_xyvtcs (GRID + OBSERVATIONS)
    # =====================================================================
    print("\n3. Building df_xyvtcs...")

    xcol, ycol, vcol = "x", "y", "value"
    tol = 1e-8

    grid_xy = df_grid_xyct[[xcol, ycol]].to_numpy()
    obs_xy = df_fl_xyvctr[[xcol, ycol]].to_numpy()
    obs_val = df_fl_xyvctr[vcol].to_numpy()

    tree = cKDTree(grid_xy)
    dist, gi = tree.query(obs_xy, distance_upper_bound=tol)
    valid = np.isfinite(dist) & (gi >= 0) & (gi < len(grid_xy))

    # Start from grid, add observation info
    df_xyvtcs = df_grid_xyct.copy()
    df_xyvtcs["set"] = np.zeros(len(df_xyvtcs), dtype=np.int8)
    df_xyvtcs[vcol] = np.nan

    # Handle duplicates: average obs that hit the same grid node
    matched = pd.DataFrame({"grid_idx": gi[valid].astype(int), vcol: obs_val[valid]})
    agg = matched.groupby("grid_idx", as_index=False)[vcol].mean()

    # Assign values and flags
    df_xyvtcs.loc[agg["grid_idx"].values, vcol] = agg[vcol].values
    df_xyvtcs.loc[agg["grid_idx"].values, "set"] = 1

    # Mark excluded nodes: cluster < 0 -> set = -1
    df_xyvtcs["cluster"] = pd.to_numeric(df_xyvtcs["cluster"], errors="coerce")
    neg_mask = df_xyvtcs["cluster"] < 0
    df_xyvtcs.loc[neg_mask & (df_xyvtcs["set"] != 1), "set"] = -1

    # Reorder columns
    first = [xcol, ycol, vcol, "trend", "cluster", "set"]
    rest = [c for c in df_xyvtcs.columns if c not in first]
    df_xyvtcs = df_xyvtcs.loc[:, first + rest]

    print(f"  Total grid points: {len(df_xyvtcs)}")
    print(f"  Observation points (set=1): {(df_xyvtcs['set'] == 1).sum()}")
    print(f"  Simulation points (set=0): {(df_xyvtcs['set'] == 0).sum()}")
    print(f"  Excluded points (set=-1): {(df_xyvtcs['set'] == -1).sum()}")

    # =====================================================================
    # 4. FIT VARIOGRAMS
    # =====================================================================
    print(f"\n4. Fitting {VARIOGRAM_MODEL} variograms (maxlag={MAXLAG})...")

    variograms = []
    for i in range(n_regions):
        V = cluster_variogram(df_fl_xyvctr, 'residual', i, MAXLAG, N_LAGS, VARIOGRAM_MODEL)
        r, s, n = V.parameters
        variograms.append(V)
        print(f"  Cluster {i}: range={r:.2f}, sill={s:.4f}, nugget={n:.4f}")

    df_gamma = build_variogram_dataframe(variograms)

    # =====================================================================
    # 5. EXTRACT GRID COORDINATES
    # =====================================================================
    x_coords = np.sort(df_xyvtcs['x'].unique())
    y_coords = np.sort(df_xyvtcs['y'].unique())

    print(f"\n5. Grid info:")
    print(f"  x_coords: {len(x_coords)} unique ({x_coords.min():.2f} to {x_coords.max():.2f})")
    print(f"  y_coords: {len(y_coords)} unique ({y_coords.min():.2f} to {y_coords.max():.2f})")
    print(f"  Expected shape: {ROWS} x {COLS} = {ROWS * COLS}")
    print(f"  Actual points: {len(df_xyvtcs)}")

    # =====================================================================
    # 6. SAVE CONFIG
    # =====================================================================
    print(f"\n6. Saving configuration...")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Save DataFrame to CSV
    csv_path = OUTPUT_DIR / f"config_{ENSEMBLE_NAME}_data.csv"
    df_xyvtcs.to_csv(csv_path, index=False)
    print(f"  Saved data CSV: {csv_path}")

    # Build variogram list
    variograms_list = df_gamma['Variogram'].tolist()

    # Build config dict
    config = {
        'data_file': f"config_{ENSEMBLE_NAME}_data.csv",
        'variograms_list': variograms_list,
        'vario': None,
        'mg_resols': MG_RESOLS,
        'grid_shape': [ROWS, COLS],
        'x_coords': x_coords.tolist(),
        'y_coords': y_coords.tolist(),
        'use_nst': False,
        'mgsim_kwargs': {
            'xx': 'x',
            'yy': 'y',
            'zz': 'residual',
            'kk': 'cluster',
            'num_points': NUM_POINTS,
            'radius': RADIUS,
            'sgs_or_krige': 'sgs',
        },
        'metadata': {
            'ensemble_name': ENSEMBLE_NAME,
            'use_subregions': True,
            'variogram_model': VARIOGRAM_MODEL,
            'maxlag': MAXLAG,
            'n_lags': N_LAGS,
            'smoothing': SMOOTHING,
            'spacing_fl': SPACING_FL,
            'random_seed': RANDOM_SEED,
        }
    }

    json_path = OUTPUT_DIR / f"config_{ENSEMBLE_NAME}.json"
    with open(json_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"  Saved config JSON: {json_path}")

    # =====================================================================
    # SUMMARY
    # =====================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Grid shape: {ROWS} x {COLS}")
    print(f"  MG resolutions: {MG_RESOLS}")
    print(f"  Variogram model: {VARIOGRAM_MODEL}")
    print(f"  Maxlag: {MAXLAG} (was 400 in v1)")
    print(f"  Smoothing: {SMOOTHING} (was 1000 in v1)")
    print(f"  Spacing FL: {SPACING_FL} (was 200 in v1)")
    print(f"  num_points: {NUM_POINTS}, radius: {RADIUS}")
    print(f"  Config: {json_path}")
    print(f"  Data: {csv_path}")
    print()
    print("Next steps:")
    print("  1. Copy config files to Sherlock:")
    print(f"     scp {json_path} {csv_path} \\")
    print(f"       jrines@login.sherlock.stanford.edu:/oak/stanford/groups/cyaolai/JoshRines/repos/multigrid-sgsim/configs/")
    print("  2. Submit batch job:")
    print("     sbatch submit_csb_mgsim_v2.sbatch")


if __name__ == '__main__':
    main()
