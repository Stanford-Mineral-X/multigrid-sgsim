#!/usr/bin/env python
"""
Combine MGSIM Batch Results

This script combines multiple partial NetCDF files from array jobs
into a single consolidated NetCDF file.

It also computes summary statistics (mean, variance, percentiles).

Usage:
    python combine_results.py --input-dir results --output mgsim_combined.nc
"""

import argparse
import numpy as np
import xarray as xr
from pathlib import Path
from datetime import datetime
import glob


def combine_netcdf_files(input_dir: str, output_path: str, pattern: str = "mgsim_realizations_*.nc"):
    """
    Combine multiple NetCDF files into one and compute statistics.

    Parameters
    ----------
    input_dir : str
        Directory containing partial NetCDF files
    output_path : str
        Path for the combined output file
    pattern : str
        Glob pattern to match input files
    """
    input_path = Path(input_dir)
    files = sorted(input_path.glob(pattern))

    if not files:
        print(f"No files found matching pattern '{pattern}' in {input_dir}")
        return

    print(f"Found {len(files)} files to combine:")
    for f in files:
        print(f"  - {f.name}")

    # Load and concatenate all files
    print("\nLoading and concatenating...")
    datasets = []
    for f in files:
        ds = xr.open_dataset(f)
        datasets.append(ds)
        print(f"  Loaded {f.name}: {len(ds.realization)} realizations")

    # Concatenate along realization dimension
    combined = xr.concat(datasets, dim='realization')

    # Sort by realization index
    combined = combined.sortby('realization')

    print(f"\nCombined dataset: {len(combined.realization)} total realizations")

    # Compute summary statistics
    print("Computing summary statistics...")

    simulated = combined['simulated']

    # Mean and variance across realizations
    mean_field = simulated.mean(dim='realization')
    var_field = simulated.var(dim='realization')
    std_field = simulated.std(dim='realization')

    # Percentiles
    p05 = simulated.quantile(0.05, dim='realization')
    p25 = simulated.quantile(0.25, dim='realization')
    p50 = simulated.quantile(0.50, dim='realization')  # median
    p75 = simulated.quantile(0.75, dim='realization')
    p95 = simulated.quantile(0.95, dim='realization')

    # Add statistics to dataset
    combined['mean'] = mean_field
    combined['variance'] = var_field
    combined['std'] = std_field
    combined['p05'] = p05.drop_vars('quantile')
    combined['p25'] = p25.drop_vars('quantile')
    combined['median'] = p50.drop_vars('quantile')
    combined['p75'] = p75.drop_vars('quantile')
    combined['p95'] = p95.drop_vars('quantile')

    # Interquartile range
    combined['iqr'] = combined['p75'] - combined['p25']

    # Update attributes
    combined.attrs['n_realizations'] = len(combined.realization)
    combined.attrs['combined_from'] = [f.name for f in files]
    combined.attrs['combined_at'] = datetime.now().isoformat()

    # Variable descriptions
    combined['simulated'].attrs['description'] = 'Simulated values from MGSIM'
    combined['mean'].attrs['description'] = 'Mean across all realizations'
    combined['variance'].attrs['description'] = 'Variance across all realizations'
    combined['std'].attrs['description'] = 'Standard deviation across all realizations'
    combined['median'].attrs['description'] = 'Median (50th percentile) across realizations'
    combined['p05'].attrs['description'] = '5th percentile across realizations'
    combined['p25'].attrs['description'] = '25th percentile across realizations'
    combined['p75'].attrs['description'] = '75th percentile across realizations'
    combined['p95'].attrs['description'] = '95th percentile across realizations'
    combined['iqr'].attrs['description'] = 'Interquartile range (p75 - p25)'

    # Save combined dataset
    print(f"\nSaving to {output_path}...")

    # Use compression
    encoding = {var: {'zlib': True, 'complevel': 4} for var in combined.data_vars}
    combined.to_netcdf(output_path, encoding=encoding)

    print("Done!")

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total realizations: {len(combined.realization)}")
    print(f"Grid shape: {combined.dims['y']} x {combined.dims['x']}")
    print(f"\nGlobal statistics:")
    print(f"  Mean of means: {combined['mean'].mean().values:.4f}")
    print(f"  Mean of stds:  {combined['std'].mean().values:.4f}")
    print(f"  Min value:     {combined['simulated'].min().values:.4f}")
    print(f"  Max value:     {combined['simulated'].max().values:.4f}")
    print(f"\nOutput file: {output_path}")
    print(f"File size: {Path(output_path).stat().st_size / 1e6:.1f} MB")

    return combined


def main():
    parser = argparse.ArgumentParser(description='Combine MGSIM batch results')
    parser.add_argument('--input-dir', type=str, default='results',
                        help='Directory containing partial NetCDF files')
    parser.add_argument('--output', type=str, default='mgsim_combined.nc',
                        help='Output combined NetCDF file path')
    parser.add_argument('--pattern', type=str, default='mgsim_realizations_*.nc',
                        help='Glob pattern to match input files')

    args = parser.parse_args()

    combine_netcdf_files(args.input_dir, args.output, args.pattern)


if __name__ == '__main__':
    main()
