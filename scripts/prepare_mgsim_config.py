#!/usr/bin/env python
"""
Prepare MGSIM Configuration for Batch Processing

This script prepares the input data and saves it to a pickle file
that can be used by run_mgsim_batch.py on HPC.

Usage:
    python prepare_mgsim_config.py

Edit the parameters below to match your dataset.
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


def prepare_synthetic_data_config():
    """
    Prepare configuration for the synthetic dataset.
    Modify this function for your specific dataset.
    """

    # =========================================================================
    # USER PARAMETERS - MODIFY THESE FOR YOUR DATASET
    # =========================================================================

    # Data paths
    data_dir = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data')
    gt_path = data_dir / 'gt_xyvc.csv'
    fl_path = data_dir / 'fl_xyvc.csv'

    # Grid parameters (should match your ground truth grid)
    rows = 100
    cols = 100

    # Multigrid resolutions (coarse to fine)
    mg_resols = [16, 8, 4, 2, 1]

    # MGSIM parameters
    mgsim_kwargs = {
        'xx': 'x',
        'yy': 'y',
        'zz': 'residual',  # or 'Nresidual' if using NST
        'kk': 'cluster',
        'num_points': 10,
        'radius': 400,
        'sgs_or_krige': 'sgs',
    }

    # NST parameters
    use_nst = False  # Set to True to use Normal Score Transform
    clip_nst = True
    clip_percentile = 99.0

    # Variogram parameters
    maxlag = 40
    n_lags = 20
    variogram_model = 'gaussian'

    # =========================================================================
    # LOAD AND PREPARE DATA
    # =========================================================================

    print("Loading data...")

    # Load ground truth (for grid coordinates)
    df_gt = pd.read_csv(gt_path)
    x_coords = np.sort(df_gt['x'].unique())
    y_coords = np.sort(df_gt['y'].unique())

    # Load flight line data (observations)
    df_fl = pd.read_csv(fl_path)
    df_fl.columns = ['x', 'y', 'value', 'cluster']
    df_fl['cluster'] = df_fl['cluster'].astype(int)

    print(f"Ground truth: {len(df_gt)} points")
    print(f"Flight lines: {len(df_fl)} observations")
    print(f"Grid shape: {rows} x {cols}")
    print(f"Clusters: {df_fl['cluster'].nunique()}")

    # =========================================================================
    # COMPUTE TREND
    # =========================================================================

    print("\nComputing trend...")

    # Simple example: use cluster means as trend
    # Replace this with your actual trend computation
    cluster_means = df_fl.groupby('cluster')['value'].mean()

    # Create full grid DataFrame
    df_grid = df_gt[['x', 'y', 'val', 'clust']].copy()
    df_grid.columns = ['x', 'y', 'value', 'cluster']
    df_grid['cluster'] = df_grid['cluster'].astype(int)

    # Assign trend based on cluster
    df_grid['trend'] = df_grid['cluster'].map(cluster_means)
    df_grid['set'] = 0  # grid points

    # Add observation points
    df_obs = df_fl.copy()
    df_obs['trend'] = df_obs['cluster'].map(cluster_means)
    df_obs['set'] = 1  # observation points

    # Combine
    df_xyvtcs = pd.concat([df_grid, df_obs], ignore_index=True)
    df_xyvtcs['residual'] = df_xyvtcs['value'] - df_xyvtcs['trend']

    print(f"Combined DataFrame: {len(df_xyvtcs)} points")
    print(f"  Grid points: {(df_xyvtcs['set'] == 0).sum()}")
    print(f"  Observation points: {(df_xyvtcs['set'] == 1).sum()}")

    # =========================================================================
    # COMPUTE VARIOGRAMS
    # =========================================================================

    print("\nComputing variograms...")

    # Prepare observation data for variogram fitting
    df_fl_var = df_fl.copy()
    df_fl_var['residual'] = df_fl_var['value'] - df_fl_var['cluster'].map(cluster_means)

    n_clusters = df_fl_var['cluster'].nunique()
    variograms = []

    for i in range(n_clusters):
        V = cluster_variogram(df_fl_var, 'residual', i, maxlag, n_lags, variogram_model)
        variograms.append(V)
        print(f"  Cluster {i}: range={V.parameters[0]:.1f}, sill={V.parameters[1]:.3f}, nugget={V.parameters[2]:.3f}")

    df_gamma = build_variogram_dataframe(variograms)
    print(f"\nVariogram DataFrame:\n{df_gamma}")

    # =========================================================================
    # PREPARE NST TRANSFORMER (if using)
    # =========================================================================

    nst_trans = None
    if use_nst:
        print("\nFitting NST transformer...")
        obs_residuals = df_xyvtcs.loc[df_xyvtcs['set'] == 1, 'residual'].values.reshape(-1, 1)
        nst_trans = QuantileTransformer(output_distribution='normal', n_quantiles=min(1000, len(obs_residuals)))
        nst_trans.fit(obs_residuals)

        # Update mgsim_kwargs for NST
        mgsim_kwargs['zz'] = 'Nresidual'
        mgsim_kwargs['clip_nst'] = clip_nst
        mgsim_kwargs['clip_percentile'] = clip_percentile

    # =========================================================================
    # BUILD CONFIG DICTIONARY
    # =========================================================================

    config = {
        'df_xyvtcs': df_xyvtcs,
        'df_gamma': df_gamma,
        'mg_resols': mg_resols,
        'grid_shape': (rows, cols),
        'x_coords': x_coords,
        'y_coords': y_coords,
        'use_nst': use_nst,
        'nst_trans': nst_trans,
        'mgsim_kwargs': mgsim_kwargs,
        # Metadata
        'metadata': {
            'data_source': str(fl_path),
            'n_observations': len(df_fl),
            'n_clusters': n_clusters,
            'variogram_model': variogram_model,
            'maxlag': maxlag,
            'n_lags': n_lags,
        }
    }

    return config


def main():
    # Prepare configuration
    config = prepare_synthetic_data_config()

    # Save to pickle
    output_path = Path(__file__).parent / 'mgsim_config.pkl'
    print(f"\nSaving configuration to {output_path}")

    with open(output_path, 'wb') as f:
        pickle.dump(config, f)

    print("Done!")

    # Print summary
    print("\n" + "=" * 60)
    print("CONFIGURATION SUMMARY")
    print("=" * 60)
    print(f"Grid shape: {config['grid_shape']}")
    print(f"Multigrid resolutions: {config['mg_resols']}")
    print(f"Use NST: {config['use_nst']}")
    print(f"MGSIM kwargs: {config['mgsim_kwargs']}")
    print(f"\nConfig saved to: {output_path}")
    print("\nTo run batch realizations:")
    print(f"  python run_mgsim_batch.py --config {output_path} --start 0 --end 100 --output results.nc")


if __name__ == '__main__':
    main()
