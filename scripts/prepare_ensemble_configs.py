#!/usr/bin/env python
"""
Prepare MGSIM Configurations for 5 Ensemble Comparison

Creates 5 configuration files for comparing:
1. Isotropic + Subregions + Dense flight lines
2. Isotropic + Global + Dense flight lines
3. Anisotropic + Subregions + Dense flight lines
4. Anisotropic + Global + Dense flight lines
5. Isotropic + Subregions + Sparse flight lines (compare with #1)

Usage:
    python prepare_ensemble_configs.py
"""

import pickle
import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src directory to path
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

from sklearn.preprocessing import QuantileTransformer
from variograms import cluster_variogram, build_variogram_dataframe
from skgstat import Variogram


# =============================================================================
# USER PARAMETERS - MODIFY THESE FOR YOUR DATASET
# =============================================================================

# Data paths
DATA_DIR = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data')
GT_PATH = DATA_DIR / 'gt_xyvc.csv'

# Flight line paths (dense and sparse for comparison)
FL_DENSE_PATH = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/multigrid-sgsim/demos/data/fl_xyvc_dense.csv')
FL_SPARSE_PATH = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/multigrid-sgsim/demos/data/fl_xyvc_sparse.csv')

# Grid parameters
ROWS = 100
COLS = 100

# Multigrid resolutions (coarse to fine)
MG_RESOLS = [16, 8, 4, 2, 1]

# MGSIM parameters
NUM_POINTS = 10
RADIUS = 400

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

# Output directory
OUTPUT_DIR = script_dir


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
    Compute trend surface.

    Parameters
    ----------
    df_fl : DataFrame
        Flight line observations
    df_gt : DataFrame
        Ground truth grid
    use_subregions : bool
        If True, use cluster-specific means
        If False, use global mean

    Returns
    -------
    df_xyvtcs : DataFrame
        Combined grid + observations with trend column
    """
    if use_subregions:
        # Cluster-specific means
        cluster_means = df_fl.groupby('cluster')['value'].mean()
    else:
        # Global mean
        global_mean = df_fl['value'].mean()
        cluster_means = pd.Series({c: global_mean for c in df_fl['cluster'].unique()})

    # Create full grid DataFrame
    df_grid = df_gt[['x', 'y', 'val', 'clust']].copy()
    df_grid.columns = ['x', 'y', 'value', 'cluster']
    df_grid['cluster'] = df_grid['cluster'].astype(int)

    if use_subregions:
        df_grid['trend'] = df_grid['cluster'].map(cluster_means)
    else:
        df_grid['trend'] = global_mean
    df_grid['set'] = 0  # grid points

    # Add observation points
    df_obs = df_fl.copy()
    if use_subregions:
        df_obs['trend'] = df_obs['cluster'].map(cluster_means)
    else:
        df_obs['trend'] = global_mean
    df_obs['set'] = 1  # observation points

    # For global case, set all clusters to 0
    if not use_subregions:
        df_grid['cluster'] = 0
        df_obs['cluster'] = 0

    # Combine
    df_xyvtcs = pd.concat([df_grid, df_obs], ignore_index=True)
    df_xyvtcs['residual'] = df_xyvtcs['value'] - df_xyvtcs['trend']

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


def create_config(df_xyvtcs, x_coords, y_coords, ensemble_name, use_subregions,
                  df_gamma=None, vario=None):
    """Create configuration dictionary.

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
    """
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

    # Fit NST transformer if needed
    nst_trans = None
    if USE_NST:
        obs_residuals = df_xyvtcs.loc[df_xyvtcs['set'] == 1, 'residual'].values.reshape(-1, 1)
        nst_trans = QuantileTransformer(output_distribution='normal', n_quantiles=min(1000, len(obs_residuals)))
        nst_trans.fit(obs_residuals)

    config = {
        'df_xyvtcs': df_xyvtcs,
        'df_gamma': df_gamma,  # For subregions (cluster_sgs)
        'vario': vario,        # For global (okrige_sgs)
        'mg_resols': MG_RESOLS,
        'grid_shape': (ROWS, COLS),
        'x_coords': x_coords,
        'y_coords': y_coords,
        'use_nst': USE_NST,
        'nst_trans': nst_trans,
        'mgsim_kwargs': mgsim_kwargs,
        'metadata': {
            'ensemble_name': ensemble_name,
            'use_subregions': use_subregions,
            'variogram_model': VARIOGRAM_MODEL,
            'maxlag': MAXLAG,
            'n_lags': N_LAGS,
        }
    }

    return config


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("MGSIM ENSEMBLE CONFIGURATION GENERATOR")
    print("=" * 70)

    # Define 5 ensembles:
    # 4 combinations of iso/aniso × subregions/global (dense flight lines)
    # + 1 with sparse flight lines for comparison
    ensembles = [
        {'name': 'iso_subregions_dense', 'use_subregions': True, 'anisotropic': False, 'fl_path': FL_DENSE_PATH},
        {'name': 'iso_global_dense', 'use_subregions': False, 'anisotropic': False, 'fl_path': FL_DENSE_PATH},
        {'name': 'aniso_subregions_dense', 'use_subregions': True, 'anisotropic': True, 'fl_path': FL_DENSE_PATH},
        {'name': 'aniso_global_dense', 'use_subregions': False, 'anisotropic': True, 'fl_path': FL_DENSE_PATH},
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

        # Create config
        config = create_config(
            df_xyvtcs, x_coords, y_coords,
            ensemble_name=ens['name'],
            use_subregions=ens['use_subregions'],
            df_gamma=df_gamma,
            vario=vario
        )

        # Save config
        output_path = OUTPUT_DIR / f"config_{ens['name']}.pkl"
        print(f"\nSaving config to {output_path}")
        with open(output_path, 'wb') as f:
            pickle.dump(config, f)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY - 5 ENSEMBLE CONFIGS CREATED")
    print("=" * 70)
    print("""
    1. config_iso_subregions_dense.pkl
       - Isotropic variograms, cluster-specific
       - Dense flight lines (spacing=3)

    2. config_iso_global_dense.pkl
       - Isotropic variograms, single global (okrige_sgs)
       - Dense flight lines

    3. config_aniso_subregions_dense.pkl
       - Anisotropic variograms, cluster-specific
       - Dense flight lines

    4. config_aniso_global_dense.pkl
       - Anisotropic variograms, single global (okrige_sgs)
       - Dense flight lines

    5. config_iso_subregions_sparse.pkl
       - Isotropic variograms, cluster-specific
       - Sparse flight lines (spacing=6)
       - Compare with #1 to assess flight line density effect
    """)

    print("To run all ensembles, use:")
    print("  sbatch submit_all_ensembles.sbatch")


if __name__ == '__main__':
    main()
