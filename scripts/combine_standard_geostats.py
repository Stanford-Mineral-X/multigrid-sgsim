#!/usr/bin/env python
"""
Combine Standard Geostatistics Results

This script combines Kriging and SGSIM results from array jobs
and computes summary statistics.

Usage:
    python combine_standard_geostats.py --input-dir results/standard_geostats/dense --output dense_combined.nc
"""

import argparse
import numpy as np
import xarray as xr
from pathlib import Path
from datetime import datetime


def combine_sgsim_files(input_dir: Path, pattern: str = "sgsim_realizations_*.nc") -> xr.Dataset:
    """Combine partial SGSIM files."""
    files = sorted(input_dir.glob(pattern))

    if not files:
        print(f"No SGSIM files found matching {pattern}")
        return None

    print(f"Found {len(files)} SGSIM files")

    datasets = []
    for f in files:
        ds = xr.open_dataset(f)
        datasets.append(ds)
        print(f"  Loaded {f.name}: {len(ds.realization)} realizations")

    # Concatenate
    combined = xr.concat(datasets, dim='realization')
    combined = combined.sortby('realization')

    # Compute statistics
    simulated = combined['simulated']
    combined['mean'] = simulated.mean(dim='realization')
    combined['variance'] = simulated.var(dim='realization')
    combined['std'] = simulated.std(dim='realization')
    combined['median'] = simulated.median(dim='realization')

    combined.attrs['n_realizations'] = len(combined.realization)
    combined.attrs['combined_at'] = datetime.now().isoformat()

    return combined


def main():
    parser = argparse.ArgumentParser(description='Combine standard geostats results')
    parser.add_argument('--input-dir', type=str, required=True,
                        help='Directory containing kriging and SGSIM results')
    parser.add_argument('--output', type=str, default='standard_geostats_combined.nc',
                        help='Output combined NetCDF file')

    args = parser.parse_args()

    input_dir = Path(args.input_dir)

    # Load kriging result
    kriging_path = input_dir / 'kriging_result.nc'
    if kriging_path.exists():
        print(f"Loading kriging result from {kriging_path}")
        ds_krige = xr.open_dataset(kriging_path)
    else:
        print(f"No kriging result found at {kriging_path}")
        ds_krige = None

    # Combine SGSIM files
    print("\nCombining SGSIM files...")
    ds_sgsim = combine_sgsim_files(input_dir)

    # Merge kriging and SGSIM
    if ds_krige is not None and ds_sgsim is not None:
        # Add kriging to SGSIM dataset
        ds_combined = ds_sgsim.copy()
        ds_combined['kriging_pred'] = ds_krige['kriging_pred']
        ds_combined['kriging_var'] = ds_krige['kriging_var']
        ds_combined.attrs['has_kriging'] = 1
    elif ds_krige is not None:
        ds_combined = ds_krige
        ds_combined.attrs['has_kriging'] = 1
    elif ds_sgsim is not None:
        ds_combined = ds_sgsim
        ds_combined.attrs['has_kriging'] = 0
    else:
        print("No data found!")
        return

    # Save
    print(f"\nSaving combined results to {args.output}")
    encoding = {var: {'zlib': True, 'complevel': 4} for var in ds_combined.data_vars}
    ds_combined.to_netcdf(args.output, encoding=encoding)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    if 'realization' in ds_combined.dims:
        print(f"SGSIM realizations: {len(ds_combined.realization)}")
        print(f"  Mean of means: {ds_combined['mean'].mean().values:.4f}")
        print(f"  Mean of stds:  {ds_combined['std'].mean().values:.4f}")
    if 'kriging_pred' in ds_combined:
        print(f"Kriging: mean={ds_combined['kriging_pred'].mean().values:.4f}")
    print(f"Grid shape: {ds_combined.dims['y']} x {ds_combined.dims['x']}")


if __name__ == '__main__':
    main()
