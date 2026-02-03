#!/usr/bin/env python
"""
Run Standard Geostatistical Interpolation Methods for Comparison

This script runs Kriging and SGSIM using GStatSim for comparison with MGSIM.
Designed for use on HPC (Sherlock) with SLURM.

Methods:
- Ordinary Kriging (deterministic)
- SGSIM (stochastic, multiple realizations)

Usage:
    # Run kriging only
    python run_standard_geostats.py --config config.json --method kriging --output kriging_result.nc

    # Run SGSIM realizations
    python run_standard_geostats.py --config config.json --method sgsim --start 0 --end 100 --output sgsim_results.nc

    # Run both
    python run_standard_geostats.py --config config.json --method both --start 0 --end 100 --output-dir results/
"""

import argparse
import json
import numpy as np
import pandas as pd
import xarray as xr
import sys
import os
from pathlib import Path
from datetime import datetime
from sklearn.preprocessing import QuantileTransformer

# Add src directory to path for local imports
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

try:
    import gstatsim as gs
except ImportError:
    print("GStatSim not found. Install with: pip install gstatsim")
    sys.exit(1)

try:
    from skgstat import Variogram
except ImportError:
    print("scikit-gstat not found. Install with: pip install scikit-gstat")
    sys.exit(1)


# =============================================================================
# CONFIGURATION
# =============================================================================

def load_config(config_path: str) -> dict:
    """Load configuration from JSON file."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config


def prepare_data(config: dict) -> tuple:
    """
    Prepare data for geostatistical interpolation.

    Returns
    -------
    df_obs : DataFrame
        Observation data with columns [x, y, value, nvalue]
    pred_grid : ndarray
        Prediction grid coordinates (n_points, 2)
    grid_shape : tuple
        (rows, cols)
    x_coords, y_coords : ndarray
        1D coordinate arrays
    nst_trans : QuantileTransformer
        Fitted normal score transformer
    """
    config_dir = Path(config.get('config_dir', '.'))

    # Load data CSV
    data_file = config_dir / config['data_file']
    df = pd.read_csv(data_file)

    # Extract grid info
    grid_shape = tuple(config['grid_shape'])
    rows, cols = grid_shape
    x_coords = np.array(config['x_coords'])
    y_coords = np.array(config['y_coords'])

    # Build prediction grid (row-major order)
    xx, yy = np.meshgrid(x_coords, y_coords)
    pred_grid = np.column_stack([xx.ravel(), yy.ravel()])

    # Extract observations (set == 1)
    df_obs = df[df['set'] == 1].copy()
    df_obs = df_obs.rename(columns={'value': 'val'})

    # Fit Normal Score Transform on observations
    obs_values = df_obs['val'].values.reshape(-1, 1)
    nst_trans = QuantileTransformer(
        output_distribution='normal',
        n_quantiles=min(1000, len(obs_values))
    )
    nst_trans.fit(obs_values)

    # Transform observations
    df_obs['nval'] = nst_trans.transform(obs_values).ravel()

    print(f"Loaded {len(df_obs)} observations")
    print(f"Grid shape: {rows} x {cols} = {rows * cols} points")
    print(f"NST fitted on {len(obs_values)} values")

    return df_obs, pred_grid, grid_shape, x_coords, y_coords, nst_trans


def fit_variogram(df_obs: pd.DataFrame, config: dict) -> list:
    """
    Fit variogram to observation data.

    Returns
    -------
    vario : list
        [azimuth, nugget, major_range, minor_range, sill, vtype]
    """
    # Get variogram parameters from config or fit new
    if config.get('vario') is not None:
        # Use pre-specified variogram
        vario = config['vario']
        print(f"Using pre-specified variogram: {vario}")
        return vario

    # Fit new variogram
    maxlag = config.get('maxlag', 40)
    n_lags = config.get('n_lags', 20)
    model = config.get('variogram_model', 'gaussian')

    print(f"Fitting {model} variogram (maxlag={maxlag}, n_lags={n_lags})...")

    coords = df_obs[['x', 'y']].values
    values = df_obs['nval'].values

    V = Variogram(coords, values, maxlag=maxlag, n_lags=n_lags, model=model)

    # Extract parameters: skgstat returns [range, sill, nugget]
    range_param, sill, nugget = V.parameters

    # Build vario list for GStatSim
    # Format: [azimuth, nugget, major_range, minor_range, sill, vtype]
    azimuth = 0.0  # Isotropic
    vtype = model.capitalize()

    vario = [azimuth, nugget, range_param, range_param, sill, vtype]

    print(f"Fitted variogram: range={range_param:.2f}, sill={sill:.4f}, nugget={nugget:.4f}")

    return vario


# =============================================================================
# KRIGING
# =============================================================================

def run_kriging(df_obs: pd.DataFrame, pred_grid: np.ndarray, vario: list,
                grid_shape: tuple, x_coords: np.ndarray, y_coords: np.ndarray,
                nst_trans, config: dict) -> xr.Dataset:
    """
    Run Ordinary Kriging interpolation.

    Returns
    -------
    ds : xr.Dataset
        Dataset with kriging results
    """
    rows, cols = grid_shape
    k = config.get('num_points', 10)
    rad = config.get('radius', 40)

    print(f"\nRunning Ordinary Kriging...")
    print(f"  Neighbors: {k}")
    print(f"  Search radius: {rad}")
    print(f"  Variogram: {vario}")

    # Run kriging on NST-transformed values
    pred, var = gs.Interpolation.okrige(
        pred_grid, df_obs, 'x', 'y', 'nval', k, vario, rad
    )

    # Back-transform predictions
    pred_trans = nst_trans.inverse_transform(pred.reshape(-1, 1)).ravel()

    # Reshape to grid
    pred_grid_2d = pred_trans.reshape(rows, cols)
    var_grid_2d = var.reshape(rows, cols)
    pred_norm_grid = pred.reshape(rows, cols)

    # Create dataset
    ds = xr.Dataset(
        data_vars={
            'kriging_pred': (['y', 'x'], pred_grid_2d.astype(np.float32)),
            'kriging_var': (['y', 'x'], var_grid_2d.astype(np.float32)),
            'kriging_norm': (['y', 'x'], pred_norm_grid.astype(np.float32)),
        },
        coords={
            'y': y_coords,
            'x': x_coords,
        },
        attrs={
            'method': 'Ordinary Kriging',
            'description': 'Kriging interpolation via GStatSim',
            'created': datetime.now().isoformat(),
            'neighbors': k,
            'search_radius': rad,
            'variogram': str(vario),
        }
    )

    print(f"  Kriging complete: mean={np.nanmean(pred_trans):.3f}, var={np.nanvar(pred_trans):.3f}")

    return ds


# =============================================================================
# SGSIM
# =============================================================================

def run_sgsim(df_obs: pd.DataFrame, pred_grid: np.ndarray, vario: list,
              grid_shape: tuple, x_coords: np.ndarray, y_coords: np.ndarray,
              nst_trans, config: dict,
              start_idx: int, end_idx: int, seed_offset: int = 0) -> xr.Dataset:
    """
    Run Sequential Gaussian Simulation (SGSIM) realizations.

    Returns
    -------
    ds : xr.Dataset
        Dataset with SGSIM realizations
    """
    rows, cols = grid_shape
    k = config.get('num_points', 10)
    rad = config.get('radius', 40)
    n_realizations = end_idx - start_idx

    print(f"\nRunning SGSIM realizations {start_idx} to {end_idx - 1}...")
    print(f"  Neighbors: {k}")
    print(f"  Search radius: {rad}")
    print(f"  Variogram: {vario}")

    # Pre-allocate arrays
    realizations_norm = np.full((n_realizations, rows, cols), np.nan, dtype=np.float32)
    realizations_trans = np.full((n_realizations, rows, cols), np.nan, dtype=np.float32)

    for i, real_idx in enumerate(range(start_idx, end_idx)):
        print(f"  [{datetime.now().strftime('%H:%M:%S')}] Realization {real_idx} ({i+1}/{n_realizations})")

        # Set seed for reproducibility
        rng = np.random.default_rng(real_idx + seed_offset)

        try:
            # Run SGSIM
            sim = gs.Interpolation.okrige_sgs(
                pred_grid, df_obs, 'x', 'y', 'nval', k, vario, rad, seed=rng
            )

            # Back-transform
            sim_trans = nst_trans.inverse_transform(sim.reshape(-1, 1)).ravel()

            # Store
            realizations_norm[i] = sim.reshape(rows, cols)
            realizations_trans[i] = sim_trans.reshape(rows, cols)

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    # Create dataset
    realization_indices = np.arange(start_idx, end_idx)

    ds = xr.Dataset(
        data_vars={
            'simulated': (['realization', 'y', 'x'], realizations_trans),
            'simulated_norm': (['realization', 'y', 'x'], realizations_norm),
        },
        coords={
            'realization': realization_indices,
            'y': y_coords,
            'x': x_coords,
        },
        attrs={
            'method': 'SGSIM (OK-based)',
            'description': 'Sequential Gaussian Simulation via GStatSim',
            'created': datetime.now().isoformat(),
            'start_realization': start_idx,
            'end_realization': end_idx,
            'n_realizations': n_realizations,
            'neighbors': k,
            'search_radius': rad,
            'variogram': str(vario),
        }
    )

    print(f"  SGSIM complete: {n_realizations} realizations")

    return ds


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Run standard geostatistical methods')
    parser.add_argument('--config', type=str, required=True,
                        help='Path to JSON configuration file')
    parser.add_argument('--method', type=str, choices=['kriging', 'sgsim', 'both'],
                        default='both', help='Method to run')
    parser.add_argument('--start', type=int, default=0,
                        help='Starting realization index for SGSIM')
    parser.add_argument('--end', type=int, default=100,
                        help='Ending realization index for SGSIM')
    parser.add_argument('--output', type=str,
                        help='Output file path (for single method)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory (for both methods)')
    parser.add_argument('--seed-offset', type=int, default=0,
                        help='Offset for random seeds')

    args = parser.parse_args()

    # Load configuration
    print(f"Loading configuration from {args.config}")
    config = load_config(args.config)
    config['config_dir'] = str(Path(args.config).parent)

    # Prepare data
    print("\nPreparing data...")
    df_obs, pred_grid, grid_shape, x_coords, y_coords, nst_trans = prepare_data(config)

    # Fit variogram
    vario = fit_variogram(df_obs, config)

    # Run methods
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.method in ['kriging', 'both']:
        ds_krige = run_kriging(
            df_obs, pred_grid, vario, grid_shape,
            x_coords, y_coords, nst_trans, config
        )

        if args.method == 'kriging' and args.output:
            output_path = args.output
        else:
            output_path = output_dir / 'kriging_result.nc'

        print(f"\nSaving kriging results to {output_path}")
        ds_krige.to_netcdf(output_path)

    if args.method in ['sgsim', 'both']:
        ds_sgsim = run_sgsim(
            df_obs, pred_grid, vario, grid_shape,
            x_coords, y_coords, nst_trans, config,
            args.start, args.end, args.seed_offset
        )

        if args.method == 'sgsim' and args.output:
            output_path = args.output
        else:
            output_path = output_dir / f'sgsim_realizations_{args.start}_{args.end}.nc'

        print(f"\nSaving SGSIM results to {output_path}")
        ds_sgsim.to_netcdf(output_path)

    print("\nDone!")


if __name__ == '__main__':
    main()
