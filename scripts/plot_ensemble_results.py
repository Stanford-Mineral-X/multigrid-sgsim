#!/usr/bin/env python
"""
Plot Ensemble Results from Combined NetCDF Files

This script generates publication-quality figures from already-combined
MGSIM ensemble results.

Usage:
    python plot_ensemble_results.py --results-dir results/ --ground-truth gt_xyvc.csv --output-dir figures/
"""

import argparse
import sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
from pathlib import Path

# Add src directory to path for geosoft colormap
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

from utils import geosoft_cmap_k65


# Fixed color scales
VMIN, VMAX = -1000, 1500  # For mean fields
VAR_VMIN, VAR_VMAX = 0, 140000  # For variance fields


def load_ground_truth(gt_path: str) -> tuple:
    """Load ground truth and return as 2D array with coordinates."""
    df = pd.read_csv(gt_path)
    x_unique = np.sort(df['x'].unique())
    y_unique = np.sort(df['y'].unique())
    rows, cols = len(y_unique), len(x_unique)
    values = df['val'].values.reshape(rows, cols)
    extent = [x_unique.min(), x_unique.max(), y_unique.min(), y_unique.max()]
    return values, extent, x_unique, y_unique


def plot_stacked_comparison(datasets: dict, ground_truth: np.ndarray, extent: list,
                            output_path: Path, title: str = None):
    """
    Plot ground truth and ensemble means + variance stacked vertically.
    Left column: mean fields, Right column: variance fields.
    Two colorbars on the far right.

    Parameters
    ----------
    datasets : dict
        Dictionary of {name: xr.Dataset} with 'mean' and 'variance' variables
    ground_truth : np.ndarray
        Ground truth 2D array
    extent : list
        [xmin, xmax, ymin, ymax]
    output_path : Path
        Output file path
    title : str, optional
        Figure title
    """
    from matplotlib.gridspec import GridSpec

    geosoft_cmap = geosoft_cmap_k65()
    ls = LightSource(azdeg=315, altdeg=45)
    norm = Normalize(vmin=VMIN, vmax=VMAX)
    var_norm = Normalize(vmin=VAR_VMIN, vmax=VAR_VMAX)

    n_panels = 1 + len(datasets)

    # Create figure with GridSpec: mean column, variance column, two colorbars
    fig = plt.figure(figsize=(24, 4 * n_panels))
    gs = GridSpec(n_panels, 4, width_ratios=[1, 1, 0.03, 0.03], wspace=0.08)

    # Ground truth (no variance for ground truth, leave variance panel empty or show zeros)
    ax = fig.add_subplot(gs[0, 0])
    gt_rgb = ls.shade(ground_truth, cmap=geosoft_cmap, blend_mode='soft', vmin=VMIN, vmax=VMAX)
    ax.imshow(gt_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax.set_title('Ground Truth', fontsize=12)
    ax.set_xlabel('x'); ax.set_ylabel('y')

    # Empty variance panel for ground truth row
    ax_var = fig.add_subplot(gs[0, 1])
    ax_var.text(0.5, 0.5, 'N/A\n(Ground Truth)', ha='center', va='center',
                transform=ax_var.transAxes, fontsize=14, color='gray')
    ax_var.set_title('Variance', fontsize=12)
    ax_var.axis('off')

    # Ensemble means and variances
    for i, (name, ds) in enumerate(datasets.items()):
        # Mean field
        ax = fig.add_subplot(gs[i + 1, 0])
        mean_field = ds['mean'].values
        mean_rgb = ls.shade(mean_field, cmap=geosoft_cmap, blend_mode='soft', vmin=VMIN, vmax=VMAX)
        ax.imshow(mean_rgb, origin='lower', extent=extent, interpolation='nearest')
        ax.set_title(f'{name} - Mean', fontsize=12)
        ax.set_xlabel('x'); ax.set_ylabel('y')

        # Variance field
        ax_var = fig.add_subplot(gs[i + 1, 1])
        if 'variance' in ds:
            var_field = ds['variance'].values
        elif 'std' in ds:
            var_field = ds['std'].values ** 2
        else:
            var_field = np.zeros_like(mean_field)

        var_rgb = ls.shade(var_field, cmap=geosoft_cmap, blend_mode='soft', vmin=VAR_VMIN, vmax=VAR_VMAX)
        ax_var.imshow(var_rgb, origin='lower', extent=extent, interpolation='nearest')
        ax_var.set_title(f'{name} - Variance', fontsize=12)
        ax_var.set_xlabel('x'); ax_var.set_ylabel('y')

    # Colorbar for mean values
    cbar_ax1 = fig.add_subplot(gs[:, 2])
    sm1 = plt.cm.ScalarMappable(norm=norm, cmap=geosoft_cmap)
    sm1.set_array([])
    cbar1 = fig.colorbar(sm1, cax=cbar_ax1, orientation='vertical')
    cbar1.set_label('Value')

    # Colorbar for variance
    cbar_ax2 = fig.add_subplot(gs[:, 3])
    sm2 = plt.cm.ScalarMappable(norm=var_norm, cmap=geosoft_cmap)
    sm2.set_array([])
    cbar2 = fig.colorbar(sm2, cax=cbar_ax2, orientation='vertical')
    cbar2.set_label('Variance')

    if title:
        fig.suptitle(title, fontsize=14, y=1.01)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_realizations(ds: xr.Dataset, extent: list, output_path: Path,
                      name: str, n_realizations: int = 4):
    """
    Plot random realizations, then mean field, then variance field stacked vertically.
    Order: realizations on top, mean below, variance at bottom.
    Colorbars next to mean and variance fields only.
    """
    from matplotlib.gridspec import GridSpec

    geosoft_cmap = geosoft_cmap_k65()
    ls = LightSource(azdeg=315, altdeg=45)
    norm = Normalize(vmin=VMIN, vmax=VMAX)
    var_norm = Normalize(vmin=VAR_VMIN, vmax=VAR_VMAX)

    # Select random realizations
    n_real = len(ds.realization)
    np.random.seed(42)
    random_indices = np.random.choice(n_real, size=min(n_realizations, n_real), replace=False)
    random_indices.sort()

    # Rows: realizations first, then mean, then variance
    n_real_panels = len(random_indices)
    n_total_panels = n_real_panels + 2  # +2 for mean and variance

    # Create figure with GridSpec
    fig = plt.figure(figsize=(18, 4 * n_total_panels))
    # Main column for plots, narrow column for colorbars
    gs = GridSpec(n_total_panels, 2, width_ratios=[1, 0.03], wspace=0.05)

    # Get mean and variance fields
    mean_field = ds['mean'].values
    if 'variance' in ds:
        var_field = ds['variance'].values
    elif 'std' in ds:
        var_field = ds['std'].values ** 2
    else:
        var_field = np.zeros_like(mean_field)

    # Plot realizations first (no individual colorbars)
    for i, real_idx in enumerate(random_indices):
        ax = fig.add_subplot(gs[i, 0])
        realization = ds['simulated'].isel(realization=real_idx).values
        real_rgb = ls.shade(realization, cmap=geosoft_cmap, blend_mode='soft', vmin=VMIN, vmax=VMAX)
        ax.imshow(real_rgb, origin='lower', extent=extent, interpolation='nearest')
        ax.set_title(f'Realization {int(real_idx)}', fontsize=12)
        ax.set_xlabel('x'); ax.set_ylabel('y')

    # Mean field with its colorbar
    mean_row = n_real_panels
    ax_mean = fig.add_subplot(gs[mean_row, 0])
    mean_rgb = ls.shade(mean_field, cmap=geosoft_cmap, blend_mode='soft', vmin=VMIN, vmax=VMAX)
    ax_mean.imshow(mean_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax_mean.set_title(f'{name} - Ensemble Mean (n={n_real})', fontsize=12)
    ax_mean.set_xlabel('x'); ax_mean.set_ylabel('y')

    # Colorbar for mean (next to mean row only)
    cbar_ax_mean = fig.add_subplot(gs[mean_row, 1])
    sm_mean = plt.cm.ScalarMappable(norm=norm, cmap=geosoft_cmap)
    sm_mean.set_array([])
    cbar_mean = fig.colorbar(sm_mean, cax=cbar_ax_mean, orientation='vertical')
    cbar_mean.set_label('Value')

    # Variance field with its colorbar (using geosoft colormap)
    var_row = n_real_panels + 1
    ax_var = fig.add_subplot(gs[var_row, 0])
    var_rgb = ls.shade(var_field, cmap=geosoft_cmap, blend_mode='soft', vmin=VAR_VMIN, vmax=VAR_VMAX)
    ax_var.imshow(var_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax_var.set_title(f'{name} - Ensemble Variance', fontsize=12)
    ax_var.set_xlabel('x'); ax_var.set_ylabel('y')

    # Colorbar for variance (next to variance row only)
    cbar_ax_var = fig.add_subplot(gs[var_row, 1])
    sm_var = plt.cm.ScalarMappable(norm=var_norm, cmap=geosoft_cmap)
    sm_var.set_array([])
    cbar_var = fig.colorbar(sm_var, cax=cbar_ax_var, orientation='vertical')
    cbar_var.set_label('Variance')

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_pairwise_differences(datasets: dict, extent: list, output_path: Path):
    """
    Plot pixel-wise absolute differences between all pairs of ensemble means.

    Parameters
    ----------
    datasets : dict
        Dictionary of {name: xr.Dataset} with 'mean' variable
    extent : list
        [xmin, xmax, ymin, ymax]
    output_path : Path
        Output file path
    """
    from matplotlib.gridspec import GridSpec
    from itertools import combinations

    ls = LightSource(azdeg=315, altdeg=45)

    # Get all pairs
    names = list(datasets.keys())
    pairs = list(combinations(names, 2))

    if len(pairs) == 0:
        print("  Need at least 2 datasets for pairwise comparison")
        return

    n_pairs = len(pairs)

    # Compute all differences to find common scale
    differences = {}
    for name1, name2 in pairs:
        mean1 = datasets[name1]['mean'].values
        mean2 = datasets[name2]['mean'].values
        diff = np.abs(mean1 - mean2)
        differences[(name1, name2)] = diff

    # Find max difference for common scale
    all_diffs = np.concatenate([d.flatten() for d in differences.values()])
    diff_max = np.nanpercentile(all_diffs, 99)
    diff_norm = Normalize(vmin=0, vmax=diff_max)

    # Use a sequential colormap for absolute differences
    diff_cmap = plt.cm.hot_r

    # Create figure with GridSpec
    fig = plt.figure(figsize=(18, 4 * n_pairs))
    gs = GridSpec(n_pairs, 2, width_ratios=[1, 0.03], wspace=0.05)

    for i, (name1, name2) in enumerate(pairs):
        ax = fig.add_subplot(gs[i, 0])
        diff = differences[(name1, name2)]

        # Apply hillshading
        diff_rgb = ls.shade(diff, cmap=diff_cmap, blend_mode='soft', vmin=0, vmax=diff_max)
        ax.imshow(diff_rgb, origin='lower', extent=extent, interpolation='nearest')

        # Compute stats
        mean_diff = np.nanmean(diff)
        max_diff = np.nanmax(diff)

        ax.set_title(f'|{name1} - {name2}|\nMean: {mean_diff:.1f}, Max: {max_diff:.1f}', fontsize=11)
        ax.set_xlabel('x'); ax.set_ylabel('y')

    # Single colorbar spanning all rows
    cbar_ax = fig.add_subplot(gs[:, 1])
    sm = plt.cm.ScalarMappable(norm=diff_norm, cmap=diff_cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Absolute Difference')

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Plot ensemble results from combined NC files')
    parser.add_argument('--results-dir', type=str, required=True,
                        help='Directory containing combined NC files')
    parser.add_argument('--ground-truth', type=str, required=True,
                        help='Path to ground truth CSV')
    parser.add_argument('--output-dir', type=str, default='figures',
                        help='Output directory for figures')
    parser.add_argument('--pattern', type=str, default='*_combined.nc',
                        help='Glob pattern for combined NC files')

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load ground truth
    print("Loading ground truth...")
    ground_truth, extent, x_coords, y_coords = load_ground_truth(args.ground_truth)
    print(f"  Grid shape: {ground_truth.shape}")

    # Find combined NC files
    nc_files = sorted(results_dir.glob(args.pattern))
    if not nc_files:
        # Try looking in subdirectories
        nc_files = sorted(results_dir.glob(f'**/{args.pattern}'))

    if not nc_files:
        print(f"No files matching '{args.pattern}' found in {results_dir}")
        return

    print(f"\nFound {len(nc_files)} combined NC files:")
    for f in nc_files:
        print(f"  {f.name}")

    # Load datasets
    datasets = {}
    for f in nc_files:
        name = f.stem.replace('_combined', '')
        print(f"\nLoading {name}...")
        ds = xr.open_dataset(f)
        datasets[name] = ds
        if 'realization' in ds.dims:
            print(f"  {len(ds.realization)} realizations")

    # Generate comparison figure (all ensemble means + ground truth)
    print("\nGenerating comparison figure...")
    plot_stacked_comparison(
        datasets, ground_truth, extent,
        output_dir / 'ensemble_comparison.png',
        title='MGSIM Ensemble Comparison'
    )

    # Generate realization figures for each ensemble
    print("\nGenerating realization figures...")
    for name, ds in datasets.items():
        if 'realization' in ds.dims:
            plot_realizations(
                ds, extent,
                output_dir / f'{name}_realizations.png',
                name
            )

    # Generate pairwise absolute difference figure
    if len(datasets) >= 2:
        print("\nGenerating pairwise difference figure...")
        plot_pairwise_differences(
            datasets, extent,
            output_dir / 'pairwise_differences.png'
        )

    print(f"\nAll figures saved to {output_dir}")


if __name__ == '__main__':
    main()
