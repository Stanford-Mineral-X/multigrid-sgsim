#!/usr/bin/env python
"""
Prepare MGSIM Configuration for CSB with 6 flight lines held out.

Matches demo_csb_matern_cluster_nst.ipynb:
  - NST on observation residuals (QuantileTransformer)
  - Per-cluster variograms fitted on NST'd residuals (spherical → Matern)
  - mgsim_nst with clip_nst=True, clip_percentile=99.0

Usage:
    python generate_csb_6linesout.py                  # first, create the holdout dataset
    python prepare_csb_matern_nst_6linesout_config.py # then, create the config
"""

import json
import numpy as np
import pandas as pd
import sys
from pathlib import Path
from scipy.spatial import cKDTree
from sklearn.preprocessing import QuantileTransformer

# Add src directory to path
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

from variograms import cluster_variogram, build_variogram_dataframe
import trendmaking


# =============================================================================
# USER PARAMETERS  (must match demo_csb_matern_cluster_nst.ipynb)
# =============================================================================

# Data path (uses 6-lines-out holdout version)
DATA_DIR = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data')
CSB_PATH = DATA_DIR / 'csb_xyvmmc_6linesout.csv'

# Grid parameters
ROWS = 248
COLS = 426

# Multigrid resolutions (notebook: [500, 250, 50])
MG_RESOLS = [500, 250, 50]

# MGSIM parameters (notebook: num_points=8, radius=600)
NUM_POINTS = 8
RADIUS = 600

# Trend parameters (notebook: smoothing=10000, spacing_fl=1000)
SMOOTHING = 10000.0
SPACING_FL = 1000

# Variogram fitting parameters (fit spherical, convert to Matern)
MAXLAG = 2000
N_LAGS = 30
FIT_MODEL = 'spherical'

# Matern smoothness (notebook cell-3: matern_s = 1.0)
MATERN_S = 1.0

# NST parameters
NST_N_QUANTILES = 1000
CLIP_NST = True
CLIP_PERCENTILE = 99.0

# Output
OUTPUT_DIR = script_dir.parent / 'configs'
ENSEMBLE_NAME = 'csb_mat_nst_6linesout'

# Random seed
RANDOM_SEED = 42


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("CSB MGSIM CONFIG — MATERN + NST — 6 FLIGHT LINES HELD OUT")
    print(f"  Matern s={MATERN_S}, NST clip={CLIP_PERCENTILE}%")
    print(f"  mg_resols={MG_RESOLS}, num_points={NUM_POINTS}, radius={RADIUS}")
    print("=" * 70)

    np.random.seed(RANDOM_SEED)

    # =====================================================================
    # 1. LOAD DATA
    # =====================================================================
    print("\n1. Loading CSB holdout data...")
    df_xyvmmc = pd.read_csv(CSB_PATH)

    # Extract retained flight lines (fl_mask == 1)
    df_fl_only = df_xyvmmc[df_xyvmmc["fl_mask"] == 1].drop(
        columns=["fl_mask", "grid_mask"], errors="ignore"
    )
    fl_xyvc = df_fl_only.values

    # Grid (all points, same as original)
    df_grid_xyvc = df_xyvmmc.drop(columns=["fl_mask", "grid_mask"])
    grid_xyvc = df_grid_xyvc.values
    grid_xyc = np.column_stack([grid_xyvc[:, :2], grid_xyvc[:, 3]])

    n_regions = len(np.unique(fl_xyvc[:, 3]))
    print(f"  Flight line points (retained): {len(fl_xyvc)}")
    print(f"  Grid points: {len(grid_xyvc)}")
    print(f"  Clusters: {n_regions}")

    # =====================================================================
    # 2. COMPUTE TREND
    # =====================================================================
    print("\n2. Computing RBF trend...")
    fl_xyvct, grid_xyct = trendmaking.make_trend(
        fl_xyvc, grid_xyc, SMOOTHING, SPACING_FL
    )

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

    print(f"  Residual mean: {df_fl_xyvctr['residual'].mean():.1f}")
    print(f"  Residual std:  {df_fl_xyvctr['residual'].std():.1f}")

    # =====================================================================
    # 3. FIT NST ON OBSERVATION RESIDUALS
    # =====================================================================
    print("\n3. Fitting NST on observation residuals...")

    obs_residuals = df_fl_xyvctr['residual'].values.reshape(-1, 1)
    n_quantiles = min(NST_N_QUANTILES, len(obs_residuals))

    nst_trans = QuantileTransformer(
        output_distribution='normal',
        n_quantiles=n_quantiles
    )
    nst_trans.fit(obs_residuals)

    nst_residuals = nst_trans.transform(obs_residuals).ravel()
    df_fl_xyvctr['nst_residual'] = nst_residuals

    print(f"  NST fitted on {len(obs_residuals)} observations (n_quantiles={n_quantiles})")
    print(f"  Raw residuals:  mean={obs_residuals.mean():.1f}, std={obs_residuals.std():.1f}")
    print(f"  NST residuals:  mean={nst_residuals.mean():.3f}, std={nst_residuals.std():.3f}")
    print(f"  NST range: [{nst_residuals.min():.2f}, {nst_residuals.max():.2f}]")

    # =====================================================================
    # 4. FIT PER-CLUSTER VARIOGRAMS ON NST'd RESIDUALS
    # =====================================================================
    print(f"\n4. Fitting {FIT_MODEL} variograms on NST residuals (maxlag={MAXLAG})...")

    clusters = sorted(df_fl_xyvctr['cluster'].unique().astype(int))
    cluster_vgrams = []
    for cl in clusters:
        V = cluster_variogram(
            df_fl_xyvctr, 'nst_residual', cl,
            maxlag=MAXLAG, n_lags=N_LAGS, model=FIT_MODEL
        )
        r, s, n = V.parameters
        n_pts = (df_fl_xyvctr['cluster'] == cl).sum()
        print(f"  Cluster {cl:2d} ({n_pts:5d} pts): range={r:.0f}, sill={s:.4f}, nugget={n:.4f}")
        cluster_vgrams.append(V)

    # Build df_gamma with spherical fit, then convert to Matern
    df_gamma_sph = build_variogram_dataframe(cluster_vgrams, vtype='Spherical')

    df_gamma = df_gamma_sph.copy()
    df_gamma['Variogram'] = df_gamma['Variogram'].apply(
        lambda v: [v[0], v[1], v[2], v[3], v[4], 'Matern', MATERN_S]
    )

    print(f"\n  Converted to Matern (s={MATERN_S})")

    # =====================================================================
    # 5. BUILD df_xyvtcs (GRID + OBSERVATIONS)
    # =====================================================================
    print("\n5. Building df_xyvtcs...")

    xcol, ycol, vcol = "x", "y", "value"
    tol = 1e-8

    grid_xy = df_grid_xyct[[xcol, ycol]].to_numpy()
    obs_xy = df_fl_xyvctr[[xcol, ycol]].to_numpy()
    obs_val = df_fl_xyvctr[vcol].to_numpy()

    tree = cKDTree(grid_xy)
    dist, gi = tree.query(obs_xy, distance_upper_bound=tol)
    valid = np.isfinite(dist) & (gi >= 0) & (gi < len(grid_xy))

    df_xyvtcs = df_grid_xyct.copy()
    df_xyvtcs["set"] = np.zeros(len(df_xyvtcs), dtype=np.int8)
    df_xyvtcs[vcol] = np.nan

    matched = pd.DataFrame({"grid_idx": gi[valid].astype(int), vcol: obs_val[valid]})
    agg = matched.groupby("grid_idx", as_index=False)[vcol].mean()

    df_xyvtcs.loc[agg["grid_idx"].values, vcol] = agg[vcol].values
    df_xyvtcs.loc[agg["grid_idx"].values, "set"] = 1

    df_xyvtcs["cluster"] = pd.to_numeric(df_xyvtcs["cluster"], errors="coerce")
    neg_mask = df_xyvtcs["cluster"] < 0
    df_xyvtcs.loc[neg_mask & (df_xyvtcs["set"] != 1), "set"] = -1

    first = [xcol, ycol, vcol, "trend", "cluster", "set"]
    rest = [c for c in df_xyvtcs.columns if c not in first]
    df_xyvtcs = df_xyvtcs.loc[:, first + rest]

    print(f"  Observation points (set=1): {(df_xyvtcs['set'] == 1).sum()}")
    print(f"  Simulation points (set=0): {(df_xyvtcs['set'] == 0).sum()}")
    print(f"  Excluded points (set=-1): {(df_xyvtcs['set'] == -1).sum()}")

    # =====================================================================
    # 6. EXTRACT GRID COORDINATES
    # =====================================================================
    x_coords = np.sort(df_xyvtcs['x'].unique())
    y_coords = np.sort(df_xyvtcs['y'].unique())

    # =====================================================================
    # 7. SAVE CONFIG
    # =====================================================================
    print(f"\n7. Saving configuration...")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Save DataFrame to CSV
    csv_path = OUTPUT_DIR / f"config_{ENSEMBLE_NAME}_data.csv"
    df_xyvtcs.to_csv(csv_path, index=False)
    print(f"  Saved data CSV: {csv_path}")

    # Build variogram list (Matern format: [az, nugget, range, range, sill, 'Matern', s])
    variograms_list = df_gamma['Variogram'].tolist()

    config = {
        'data_file': f"config_{ENSEMBLE_NAME}_data.csv",
        'variograms_list': variograms_list,
        'vario': None,
        'mg_resols': MG_RESOLS,
        'grid_shape': [ROWS, COLS],
        'x_coords': x_coords.tolist(),
        'y_coords': y_coords.tolist(),
        'use_nst': True,
        'nst_n_quantiles': n_quantiles,
        'mgsim_kwargs': {
            'xx': 'x',
            'yy': 'y',
            'zz': 'Nresidual',
            'kk': 'cluster',
            'num_points': NUM_POINTS,
            'radius': RADIUS,
            'sgs_or_krige': 'sgs',
            'clip_nst': CLIP_NST,
            'clip_percentile': CLIP_PERCENTILE,
        },
        'metadata': {
            'ensemble_name': ENSEMBLE_NAME,
            'use_subregions': True,
            'variogram_model': 'Matern',
            'matern_smoothness': MATERN_S,
            'fit_model': FIT_MODEL,
            'maxlag': MAXLAG,
            'n_lags': N_LAGS,
            'smoothing': SMOOTHING,
            'spacing_fl': SPACING_FL,
            'random_seed': RANDOM_SEED,
            'holdout': {
                'x_values': [548400, 551000, 553650, 555050, 558610, 560950],
                'tolerance': 55,
                'source_file': 'csb_xyvmmc_6linesout.csv',
                'validation_file': 'csb_holdout_6lines.csv',
            },
        }
    }

    json_path = OUTPUT_DIR / f"config_{ENSEMBLE_NAME}.json"
    with open(json_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"  Saved config JSON: {json_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Ensemble: {ENSEMBLE_NAME}")
    print(f"  Grid shape: {ROWS} x {COLS}")
    print(f"  MG resolutions: {MG_RESOLS}")
    print(f"  Variogram: Matern (s={MATERN_S}), fitted as {FIT_MODEL}")
    print(f"  NST: QuantileTransformer (n_quantiles={n_quantiles})")
    print(f"  Clipping: {CLIP_PERCENTILE}% percentile")
    print(f"  num_points={NUM_POINTS}, radius={RADIUS}")
    print(f"  Holdout: 6 lines at x=[548400, 551000, 553650, 555050, 558610, 560950]")
    print(f"  Config: {json_path}")
    print(f"  Data: {csv_path}")
    print()
    print("Next steps:")
    print("  1. Copy config files to Sherlock:")
    print(f"     scp {json_path} {csv_path} \\")
    print(f"       jrines@login.sherlock.stanford.edu:/oak/stanford/groups/cyaolai/JoshRines/repos/multigrid-sgsim/configs/")
    print("  2. Submit batch job:")
    print("     sbatch submit_csb_matern_nst_6linesout_mgsim.sbatch")


if __name__ == '__main__':
    main()
