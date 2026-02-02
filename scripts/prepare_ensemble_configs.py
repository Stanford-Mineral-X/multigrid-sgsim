#!/usr/bin/env python
"""
Prepare MGSIM Configurations for 6 Ensemble Comparison

Creates 6 configuration files for comparing:
1. Isotropic + Subregions + Dense flight lines (spacing=4, gap=3)
2. Isotropic + Global + Dense flight lines
3. Anisotropic + Subregions + Dense flight lines
4. Anisotropic + Global + Dense flight lines
5. Isotropic + Subregions + Medium flight lines (spacing=6, gap=5)
6. Isotropic + Subregions + Sparse flight lines (spacing=8, gap=7)

Usage:
    python prepare_ensemble_configs.py
"""

import json
import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src directory to path
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

# QuantileTransformer imported only if USE_NST is True (not currently used)
from variograms import cluster_variogram, build_variogram_dataframe
from skgstat import Variogram
import trendmaking


# =============================================================================
# USER PARAMETERS - MODIFY THESE FOR YOUR DATASET
# =============================================================================

# Data paths
DATA_DIR = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data')
GT_PATH = DATA_DIR / 'gt_xyvc.csv'

# Flight line paths (dense, medium, sparse for comparison)
FL_DENSE_PATH = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/multigrid-sgsim/demos/data/fl_xyvc_dense.csv')
FL_MEDIUM_PATH = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/multigrid-sgsim/demos/data/fl_xyvc_medium.csv')
FL_SPARSE_PATH = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/multigrid-sgsim/demos/data/fl_xyvc_sparse.csv')

# Grid parameters
ROWS = 50
COLS = 200

# Multigrid resolutions (coarse to fine)
MG_RESOLS = [16, 8, 4, 2, 1]

# MGSIM parameters
NUM_POINTS = 10
RADIUS = 400

# Trend parameters (RBF interpolation)
SMOOTHING = 100.0
LINESPACING = 5

# Variogram fitting parameters
MAXLAG = 40
N_LAGS = 20
VARIOGRAM_MODEL = 'gaussian'

# Directional variogram parameters (for anisotropic cases)
AZIMUTHS = [0, 45, 90, 135]  # directions to check
BANDWIDTH = 45  # angular tolerance

# NST parameters (set USE_NST=True if needed)
USE_NST = False
CLIP_NST = True
CLIP_PERCENTILE = 98.0

# Output directory for config files
OUTPUT_DIR = script_dir.parent / 'configs'


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_data(fl_path=None):
    """Load ground truth and flight line data.

    Parameters
    ----------
    fl_path : Path, optional
        Path to flight line CSV. If None, uses FL_DENSE_PATH.
    """
    print("Loading data...")

    # Load ground truth (for grid coordinates)
    df_gt = pd.read_csv(GT_PATH)
    x_coords = np.sort(df_gt['x'].unique())
    y_coords = np.sort(df_gt['y'].unique())

    # Load flight line data (observations)
    if fl_path is None:
        fl_path = FL_DENSE_PATH
    df_fl = pd.read_csv(fl_path)
    df_fl.columns = ['x', 'y', 'value', 'cluster']
    df_fl['cluster'] = df_fl['cluster'].astype(int)

    print(f"  Ground truth: {len(df_gt)} points")
    print(f"  Flight lines: {len(df_fl)} observations (from {fl_path.name})")
    print(f"  Grid shape: {ROWS} x {COLS}")
    print(f"  Clusters: {df_fl['cluster'].nunique()}")

    return df_gt, df_fl, x_coords, y_coords


def compute_trend(df_fl, df_gt, use_subregions=True):
    """
    Compute smooth RBF trend surface using trendmaking.make_trend().

    Creates a DataFrame with ONE row per grid point. Observation values
    are mapped to their nearest grid nodes using KDTree.

    Parameters
    ----------
    df_fl : DataFrame
        Flight line observations with columns: x, y, value, cluster
    df_gt : DataFrame
        Ground truth grid with columns: x, y, val, clust
    use_subregions : bool
        If True, retain cluster assignments for cluster_sgs
        If False, set all clusters to 0 for okrige_sgs

    Returns
    -------
    df_xyvtcs : DataFrame
        Grid DataFrame with columns: x, y, value, trend, cluster, set, residual
        - set=0 for grid points to simulate (value=NaN)
        - set=1 for observation points (value from flight lines)
    """
    from scipy.spatial import cKDTree

    # Prepare arrays for make_trend()
    # fl_xyvc: flight lines with x, y, value, cluster
    fl_xyvc = df_fl[['x', 'y', 'value', 'cluster']].values

    # grid_xyc: grid with x, y, cluster
    grid_xyc = df_gt[['x', 'y', 'clust']].values

    # Compute smooth RBF trend using trendmaking module
    print(f"  Computing RBF trend (smoothing={SMOOTHING}, linespacing={LINESPACING})...")
    fl_xyvct, grid_xyct = trendmaking.make_trend(fl_xyvc, grid_xyc, SMOOTHING, LINESPACING)

    # Create DataFrame from grid with trend
    # grid_xyct has columns: x, y, cluster, trend
    df_xyvtcs = pd.DataFrame({
        'x': grid_xyct[:, 0],
        'y': grid_xyct[:, 1],
        'cluster': grid_xyct[:, 2].astype(int),
        'trend': grid_xyct[:, 3],
    })

    # Initialize: all values NaN, all set=0 (to simulate)
    df_xyvtcs['value'] = np.nan
    df_xyvtcs['set'] = 0

    # For global case, set all clusters to 0
    if not use_subregions:
        df_xyvtcs['cluster'] = 0

    # Map observation values to nearest grid nodes using KDTree
    grid_xy = df_xyvtcs[['x', 'y']].values
    obs_xy = df_fl[['x', 'y']].values
    obs_values = df_fl['value'].values

    tree = cKDTree(grid_xy)

    # Find nearest grid node for each observation
    # Use tolerance based on grid spacing
    x_unique = np.sort(df_xyvtcs['x'].unique())
    y_unique = np.sort(df_xyvtcs['y'].unique())
    dx = np.abs(np.diff(x_unique)).min() if len(x_unique) > 1 else 1
    dy = np.abs(np.diff(y_unique)).min() if len(y_unique) > 1 else 1
    tol = np.sqrt(dx**2 + dy**2) * 0.5  # Half diagonal of grid cell

    dist, grid_idx = tree.query(obs_xy, distance_upper_bound=tol)

    # Filter valid matches (within tolerance)
    valid_mask = np.isfinite(dist)
    valid_grid_idx = grid_idx[valid_mask]
    valid_obs_values = obs_values[valid_mask]

    # Handle duplicates: if multiple observations map to same grid node, average them
    obs_df = pd.DataFrame({
        'grid_idx': valid_grid_idx,
        'value': valid_obs_values
    })
    agg = obs_df.groupby('grid_idx').agg({'value': 'mean'}).reset_index()

    # Set values and set=1 at observation locations
    df_xyvtcs.loc[agg['grid_idx'].values, 'value'] = agg['value'].values
    df_xyvtcs.loc[agg['grid_idx'].values, 'set'] = 1

    # Compute residuals (value - trend)
    df_xyvtcs['residual'] = df_xyvtcs['value'] - df_xyvtcs['trend']

    print(f"  Grid points: {len(df_xyvtcs)} (all will be simulated)")
    print(f"  Observations mapped: {len(agg)} (set=1, used for residual conditioning)")
    print(f"  Non-observation points: {(df_xyvtcs['set'] == 0).sum()} (set=0)")

    return df_xyvtcs


def fit_isotropic_variograms(df_fl, use_subregions=True):
    """
    Fit isotropic (omnidirectional) variograms.

    Returns
    -------
    If use_subregions:
        df_gamma : DataFrame - Variogram parameters per cluster for cluster_sgs
        vario : None
    Else:
        df_gamma : None
        vario : list - Single variogram [azimuth, nugget, major_range, minor_range, sill, vtype] for okrige_sgs
    """
    print("  Fitting isotropic variograms...")

    # Prepare data
    if use_subregions:
        cluster_means = df_fl.groupby('cluster')['value'].mean()
    else:
        global_mean = df_fl['value'].mean()

    df_var = df_fl.copy()
    if use_subregions:
        df_var['residual'] = df_var['value'] - df_var['cluster'].map(cluster_means)
    else:
        df_var['residual'] = df_var['value'] - global_mean
        df_var['cluster'] = 0  # single cluster

    n_clusters = df_var['cluster'].nunique()
    variograms = []

    for i in range(n_clusters):
        V = cluster_variogram(df_var, 'residual', i, MAXLAG, N_LAGS, VARIOGRAM_MODEL)
        variograms.append(V)
        print(f"    Cluster {i}: range={V.parameters[0]:.1f}, sill={V.parameters[1]:.3f}, nugget={V.parameters[2]:.3f}")

    if use_subregions:
        df_gamma = build_variogram_dataframe(variograms)
        return df_gamma, None
    else:
        # For global: return single vario list instead of df_gamma
        V = variograms[0]
        major_range, sill, nugget = V.parameters
        vtype = VARIOGRAM_MODEL.capitalize()
        vario = [0.0, nugget, major_range, major_range, sill, vtype]  # azimuth=0, isotropic
        print(f"    Global vario: {vario}")
        return None, vario


def fit_anisotropic_variograms(df_fl, use_subregions=True):
    """
    Fit anisotropic (directional) variograms.

    For each cluster, find the direction of maximum range and fit
    major/minor ranges.

    Returns
    -------
    If use_subregions:
        df_gamma : DataFrame - Variogram parameters per cluster for cluster_sgs
        vario : None
    Else:
        df_gamma : None
        vario : list - Single variogram [azimuth, nugget, major_range, minor_range, sill, vtype] for okrige_sgs
    """
    print("  Fitting anisotropic variograms...")

    # Prepare data
    if use_subregions:
        cluster_means = df_fl.groupby('cluster')['value'].mean()
    else:
        global_mean = df_fl['value'].mean()

    df_var = df_fl.copy()
    if use_subregions:
        df_var['residual'] = df_var['value'] - df_var['cluster'].map(cluster_means)
    else:
        df_var['residual'] = df_var['value'] - global_mean
        df_var['cluster'] = 0  # single cluster

    n_clusters = df_var['cluster'].nunique()

    # First fit omnidirectional to get base parameters
    variograms = []
    aniso_azimuths = []
    minor_ranges = []

    for cluster_idx in range(n_clusters):
        # Fit directional variograms to find anisotropy
        dir_ranges = []
        for az in AZIMUTHS:
            V = cluster_variogram(
                df_var, 'residual', cluster_idx,
                MAXLAG, N_LAGS, VARIOGRAM_MODEL,
                azimuth=az, bandwidth=BANDWIDTH
            )
            dir_ranges.append((az, V.parameters[0]))

        # Find max and min range directions
        max_az, max_range = max(dir_ranges, key=lambda x: x[1])
        min_az, min_range = min(dir_ranges, key=lambda x: x[1])

        # Fit omnidirectional for base parameters
        V_omni = cluster_variogram(df_var, 'residual', cluster_idx, MAXLAG, N_LAGS, VARIOGRAM_MODEL)
        variograms.append(V_omni)

        aniso_azimuths.append(max_az)
        minor_ranges.append(min_range)

        ratio = max_range / min_range if min_range > 0 else 1.0
        print(f"    Cluster {cluster_idx}: azimuth={max_az}°, major={max_range:.1f}, minor={min_range:.1f}, ratio={ratio:.2f}")

    if use_subregions:
        df_gamma = build_variogram_dataframe(
            variograms,
            azimuths=aniso_azimuths,
            minor_ranges=minor_ranges
        )
        return df_gamma, None
    else:
        # For global: return single vario list
        V = variograms[0]
        major_range, sill, nugget = V.parameters
        vtype = VARIOGRAM_MODEL.capitalize()
        vario = [aniso_azimuths[0], nugget, major_range, minor_ranges[0], sill, vtype]
        print(f"    Global vario: {vario}")
        return None, vario


def save_config(df_xyvtcs, x_coords, y_coords, ensemble_name, use_subregions,
                df_gamma=None, vario=None, output_dir=None):
    """Save configuration to portable files (CSV + JSON).

    Creates:
    - config_{ensemble_name}_data.csv: The DataFrame
    - config_{ensemble_name}.json: All other parameters

    Parameters
    ----------
    df_xyvtcs : DataFrame
        Combined grid + observations data
    x_coords, y_coords : arrays
        Grid coordinates
    ensemble_name : str
        Name of the ensemble
    use_subregions : bool
        Whether using cluster-specific variograms
    df_gamma : DataFrame, optional
        Variogram parameters for subregions mode (cluster_sgs)
    vario : list, optional
        Single variogram [azimuth, nugget, major_range, minor_range, sill, vtype] for global mode (okrige_sgs)
    output_dir : Path, optional
        Output directory
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    mgsim_kwargs = {
        'xx': 'x',
        'yy': 'y',
        'zz': 'Nresidual' if USE_NST else 'residual',
        'kk': 'cluster',
        'num_points': NUM_POINTS,
        'radius': RADIUS,
        'sgs_or_krige': 'sgs',
    }

    if USE_NST:
        mgsim_kwargs['clip_nst'] = CLIP_NST
        mgsim_kwargs['clip_percentile'] = CLIP_PERCENTILE

    # Convert df_gamma to list of variogram lists (portable format)
    # df_gamma has a 'Variogram' column where each entry is [azimuth, nugget, major, minor, sill, vtype]
    variograms_list = None
    if df_gamma is not None:
        variograms_list = df_gamma['Variogram'].tolist()

    # Save DataFrame as CSV (portable across pandas versions)
    csv_path = output_dir / f"config_{ensemble_name}_data.csv"
    df_xyvtcs.to_csv(csv_path, index=False)
    print(f"  Saved DataFrame to {csv_path}")

    # Build config dict with JSON-serializable types
    config = {
        'data_file': f"config_{ensemble_name}_data.csv",  # Relative path
        'variograms_list': variograms_list,  # List of [az, nug, maj, min, sill, vtype] per cluster
        'vario': vario,                      # Single variogram for global mode
        'mg_resols': MG_RESOLS,
        'grid_shape': [ROWS, COLS],
        'x_coords': x_coords.tolist(),
        'y_coords': y_coords.tolist(),
        'use_nst': USE_NST,
        'mgsim_kwargs': mgsim_kwargs,
        'metadata': {
            'ensemble_name': ensemble_name,
            'use_subregions': use_subregions,
            'variogram_model': VARIOGRAM_MODEL,
            'maxlag': MAXLAG,
            'n_lags': N_LAGS,
        }
    }

    # Save config as JSON
    json_path = output_dir / f"config_{ensemble_name}.json"
    with open(json_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"  Saved config to {json_path}")

    return config


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("MGSIM ENSEMBLE CONFIGURATION GENERATOR")
    print("=" * 70)

    # Define 6 ensembles:
    # 4 combinations of iso/aniso × subregions/global (dense flight lines)
    # + 2 with medium/sparse flight lines for density comparison
    ensembles = [
        {'name': 'iso_subregions_dense', 'use_subregions': True, 'anisotropic': False, 'fl_path': FL_DENSE_PATH},
        {'name': 'iso_global_dense', 'use_subregions': False, 'anisotropic': False, 'fl_path': FL_DENSE_PATH},
        {'name': 'aniso_subregions_dense', 'use_subregions': True, 'anisotropic': True, 'fl_path': FL_DENSE_PATH},
        {'name': 'aniso_global_dense', 'use_subregions': False, 'anisotropic': True, 'fl_path': FL_DENSE_PATH},
        {'name': 'iso_subregions_medium', 'use_subregions': True, 'anisotropic': False, 'fl_path': FL_MEDIUM_PATH},
        {'name': 'iso_subregions_sparse', 'use_subregions': True, 'anisotropic': False, 'fl_path': FL_SPARSE_PATH},
    ]

    for ens in ensembles:
        print(f"\n{'='*70}")
        print(f"Ensemble: {ens['name']}")
        print(f"  Subregions: {ens['use_subregions']}")
        print(f"  Anisotropic: {ens['anisotropic']}")
        print(f"  Flight lines: {ens['fl_path'].name}")
        print("=" * 70)

        # Load data for this ensemble
        df_gt, df_fl, x_coords, y_coords = load_data(fl_path=ens['fl_path'])

        # Compute trend
        print("\nComputing trend...")
        df_xyvtcs = compute_trend(df_fl, df_gt, use_subregions=ens['use_subregions'])
        print(f"  Total points: {len(df_xyvtcs)}")
        print(f"  Clusters: {df_xyvtcs['cluster'].nunique()}")

        # Fit variograms
        print("\nFitting variograms...")
        if ens['anisotropic']:
            df_gamma, vario = fit_anisotropic_variograms(df_fl, use_subregions=ens['use_subregions'])
        else:
            df_gamma, vario = fit_isotropic_variograms(df_fl, use_subregions=ens['use_subregions'])

        if df_gamma is not None:
            print(f"\nVariogram DataFrame:\n{df_gamma}")
        else:
            print(f"\nGlobal vario: {vario}")

        # Save config (CSV + JSON format for portability)
        print("\nSaving config...")
        save_config(
            df_xyvtcs, x_coords, y_coords,
            ensemble_name=ens['name'],
            use_subregions=ens['use_subregions'],
            df_gamma=df_gamma,
            vario=vario
        )

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY - 6 ENSEMBLE CONFIGS CREATED")
    print("=" * 70)
    print("""
    Each ensemble generates two files:
    - config_{name}.json: Configuration parameters (variograms, grid info)
    - config_{name}_data.csv: DataFrame with grid points and observations

    1. config_iso_subregions_dense.json
       - Isotropic variograms, cluster-specific
       - Dense flight lines (spacing=4, gap=3)

    2. config_iso_global_dense.json
       - Isotropic variograms, single global (okrige_sgs)
       - Dense flight lines

    3. config_aniso_subregions_dense.json
       - Anisotropic variograms, cluster-specific
       - Dense flight lines

    4. config_aniso_global_dense.json
       - Anisotropic variograms, single global (okrige_sgs)
       - Dense flight lines

    5. config_iso_subregions_medium.json
       - Isotropic variograms, cluster-specific
       - Medium flight lines (spacing=6, gap=5)
       - Compare with #1 and #6 to assess flight line density effect

    6. config_iso_subregions_sparse.json
       - Isotropic variograms, cluster-specific
       - Sparse flight lines (spacing=8, gap=7)
       - Compare with #1 to assess flight line density effect
    """)

    print("To run all ensembles, use:")
    print("  sbatch submit_all_ensembles.sbatch")


if __name__ == '__main__':
    main()
