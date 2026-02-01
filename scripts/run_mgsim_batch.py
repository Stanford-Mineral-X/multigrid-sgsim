#!/usr/bin/env python
"""
Batch MGSIM Realization Generator for HPC

Usage:
    python run_mgsim_batch.py --config config.json --start 0 --end 100 --output results_0_100.nc

This script runs multiple MGSIM realizations and saves them to a NetCDF file.
Designed for use with SLURM array jobs on HPC systems like Sherlock.
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

# Add src directory to path
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

from mgsim import mgsim, mgsim_nst


def load_config(config_path: str) -> dict:
    """
    Load configuration from JSON + CSV files.

    Parameters
    ----------
    config_path : str
        Path to the JSON config file (e.g., config_iso_subregions_dense.json)

    Returns
    -------
    dict
        Configuration dictionary with:
        - df_xyvtcs: pandas DataFrame loaded from CSV
        - df_gamma: reconstructed variogram DataFrame (if subregions mode)
        - vario: single variogram list (if global mode)
        - other config parameters
    """
    config_path = Path(config_path)

    # Load JSON config
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Load DataFrame from CSV
    data_file = config_path.parent / config['data_file']
    print(f"  Loading data from {data_file}")
    df_xyvtcs = pd.read_csv(data_file)

    # Reconstruct df_gamma from variograms_list if present
    df_gamma = None
    if config.get('variograms_list') is not None:
        df_gamma = pd.DataFrame({'Variogram': config['variograms_list']})

    # Convert lists back to numpy arrays
    x_coords = np.array(config['x_coords'])
    y_coords = np.array(config['y_coords'])
    grid_shape = tuple(config['grid_shape'])

    # Build the config dict expected by run_realizations
    return {
        'df_xyvtcs': df_xyvtcs,
        'df_gamma': df_gamma,
        'vario': config.get('vario'),
        'mg_resols': config['mg_resols'],
        'grid_shape': grid_shape,
        'x_coords': x_coords,
        'y_coords': y_coords,
        'use_nst': config.get('use_nst', False),
        'nst_trans': None,  # NST not supported in portable format
        'mgsim_kwargs': config.get('mgsim_kwargs', {}),
        'metadata': config.get('metadata', {}),
    }


def run_realizations(config: dict, start_idx: int, end_idx: int, seed_offset: int = 0) -> xr.Dataset:
    """
    Run multiple MGSIM realizations.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing:
        - df_xyvtcs: input DataFrame
        - df_gamma: variogram parameters
        - mg_resols: multigrid resolutions
        - grid_shape: (rows, cols) tuple
        - x_coords: 1D array of x coordinates
        - y_coords: 1D array of y coordinates
        - use_nst: bool, whether to use NST version
        - nst_trans: transformer object (if use_nst)
        - mgsim_kwargs: additional kwargs for mgsim function
    start_idx : int
        Starting realization index (inclusive)
    end_idx : int
        Ending realization index (exclusive)
    seed_offset : int
        Offset to add to random seeds for reproducibility

    Returns
    -------
    xr.Dataset
        Dataset with dimensions (realization, y, x) containing the simulated fields
    """
    # Extract config
    df_xyvtcs = config['df_xyvtcs']
    df_gamma = config.get('df_gamma', None)  # For subregions mode
    vario = config.get('vario', None)        # For global mode
    mg_resols = config['mg_resols']
    grid_shape = config['grid_shape']
    x_coords = config['x_coords']
    y_coords = config['y_coords']
    use_nst = config.get('use_nst', False)
    nst_trans = config.get('nst_trans', None)
    mgsim_kwargs = config.get('mgsim_kwargs', {})

    rows, cols = grid_shape
    n_realizations = end_idx - start_idx

    # Pre-allocate array for results
    realizations = np.full((n_realizations, rows, cols), np.nan, dtype=np.float32)

    print(f"Running {n_realizations} realizations ({start_idx} to {end_idx-1})")
    print(f"Grid shape: {rows} x {cols}")
    print(f"Multigrid resolutions: {mg_resols}")
    print(f"Using NST: {use_nst}")
    print(f"Mode: {'Global (okrige_sgs)' if vario is not None else 'Subregions (cluster_sgs)'}")
    print("=" * 60)

    for i, real_idx in enumerate(range(start_idx, end_idx)):
        # Set random seed for reproducibility
        np.random.seed(real_idx + seed_offset)

        print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Generating realization {real_idx} ({i+1}/{n_realizations})")

        try:
            # Run MGSIM - pass either df_gamma (subregions) or vario (global)
            if use_nst and nst_trans is not None:
                df_result = mgsim_nst(
                    mg_resols=mg_resols,
                    df_xyvtcs=df_xyvtcs.copy(),
                    df_gamma=df_gamma,
                    vario=vario,
                    nst_trans=nst_trans,
                    **mgsim_kwargs
                )
            else:
                df_result = mgsim(
                    mg_resols=mg_resols,
                    df_xyvtcs=df_xyvtcs.copy(),
                    df_gamma=df_gamma,
                    vario=vario,
                    **mgsim_kwargs
                )

            # Extract the newtrend column and reshape to grid
            values = df_result['newtrend'].values

            # Handle the case where we need to reshape
            # Assuming df_xyvtcs is ordered by (y, x) or we need to pivot
            if len(values) == rows * cols:
                grid = values.reshape(rows, cols)
            else:
                # Need to create grid from x, y coordinates
                grid = np.full((rows, cols), np.nan)
                x = df_result['x'].values
                y = df_result['y'].values
                newtrend = df_result['newtrend'].values

                # Map coordinates to grid indices
                x_min, x_max = x_coords.min(), x_coords.max()
                y_min, y_max = y_coords.min(), y_coords.max()
                dx = (x_max - x_min) / (cols - 1) if cols > 1 else 1
                dy = (y_max - y_min) / (rows - 1) if rows > 1 else 1

                col_idx = np.round((x - x_min) / dx).astype(int)
                row_idx = np.round((y - y_min) / dy).astype(int)

                # Clip to valid range
                col_idx = np.clip(col_idx, 0, cols - 1)
                row_idx = np.clip(row_idx, 0, rows - 1)

                grid[row_idx, col_idx] = newtrend

            realizations[i] = grid.astype(np.float32)

            print(f"  Completed: mean={np.nanmean(grid):.3f}, std={np.nanstd(grid):.3f}")

        except Exception as e:
            print(f"  ERROR in realization {real_idx}: {e}")
            # Leave as NaN
            continue

    # Create xarray Dataset
    realization_indices = np.arange(start_idx, end_idx)

    ds = xr.Dataset(
        data_vars={
            'simulated': (['realization', 'y', 'x'], realizations),
        },
        coords={
            'realization': realization_indices,
            'y': y_coords,
            'x': x_coords,
        },
        attrs={
            'description': 'MGSIM realizations',
            'created': datetime.now().isoformat(),
            'start_realization': start_idx,
            'end_realization': end_idx,
            'n_realizations': n_realizations,
            'mg_resols': str(mg_resols),
            'use_nst': int(use_nst),  # Convert bool to int for netCDF4 compatibility
        }
    )

    return ds


def main():
    parser = argparse.ArgumentParser(description='Run batch MGSIM realizations')
    parser.add_argument('--config', type=str, required=True,
                        help='Path to JSON configuration file')
    parser.add_argument('--start', type=int, required=True,
                        help='Starting realization index (inclusive)')
    parser.add_argument('--end', type=int, required=True,
                        help='Ending realization index (exclusive)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output NetCDF file path')
    parser.add_argument('--seed-offset', type=int, default=0,
                        help='Offset for random seeds (default: 0)')

    args = parser.parse_args()

    # Load configuration from JSON + CSV
    print(f"Loading configuration from {args.config}")
    config = load_config(args.config)

    # Run realizations
    ds = run_realizations(
        config=config,
        start_idx=args.start,
        end_idx=args.end,
        seed_offset=args.seed_offset
    )

    # Save to NetCDF
    print(f"\nSaving results to {args.output}")
    ds.to_netcdf(args.output)
    print("Done!")

    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Realizations: {args.start} to {args.end - 1}")
    print(f"Mean of means: {ds['simulated'].mean(dim=['y', 'x']).mean().values:.4f}")
    print(f"Mean of stds: {ds['simulated'].std(dim=['y', 'x']).mean().values:.4f}")


if __name__ == '__main__':
    main()
