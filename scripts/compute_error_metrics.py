#!/usr/bin/env python
"""
Compute Error Metrics Against Ground Truth

This script computes error maps and quantitative metrics comparing
MGSIM results to ground truth data.

Metrics computed:
- RMSE (Root Mean Square Error)
- MAE (Mean Absolute Error)
- R² Score (Coefficient of Determination)
- SSIM (Structural Similarity Index)

Usage:
    python compute_error_metrics.py --results mgsim_combined.nc --ground-truth gt_xyvc.csv --output error_analysis.nc
"""

import argparse
import sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
from pathlib import Path
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from skimage.metrics import structural_similarity as ssim

# Add src directory to path for geosoft colormap
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

from utils import geosoft_cmap_k65


def load_ground_truth(gt_path: str, rows: int, cols: int) -> np.ndarray:
    """Load ground truth and reshape to grid."""
    df = pd.read_csv(gt_path)
    values = df['val'].values
    return values.reshape(rows, cols)


def compute_metrics(pred: np.ndarray, true: np.ndarray) -> dict:
    """
    Compute error metrics between predicted and true fields.

    Parameters
    ----------
    pred : np.ndarray
        Predicted values (2D array)
    true : np.ndarray
        Ground truth values (2D array)

    Returns
    -------
    dict
        Dictionary with RMSE, MAE, R², and SSIM
    """
    # Flatten for sklearn metrics
    pred_flat = pred.flatten()
    true_flat = true.flatten()

    # Handle NaN
    mask = ~(np.isnan(pred_flat) | np.isnan(true_flat))
    pred_valid = pred_flat[mask]
    true_valid = true_flat[mask]

    # RMSE
    rmse = np.sqrt(mean_squared_error(true_valid, pred_valid))

    # MAE
    mae = mean_absolute_error(true_valid, pred_valid)

    # R² Score
    r2 = r2_score(true_valid, pred_valid)

    # Bias (mean error)
    bias = np.mean(pred_valid - true_valid)

    # SSIM (normalize to 0-1 range)
    pred_norm = (pred - np.nanmin(pred)) / (np.nanmax(pred) - np.nanmin(pred) + 1e-10)
    true_norm = (true - np.nanmin(true)) / (np.nanmax(true) - np.nanmin(true) + 1e-10)
    ssim_val = ssim(true_norm, pred_norm, data_range=1.0)

    return {
        'rmse': rmse,
        'mae': mae,
        'r2': r2,
        'bias': bias,
        'ssim': ssim_val,
    }


def analyze_realizations(results_path: str, gt_path: str, output_path: str = None,
                         figure_path: str = None):
    """
    Analyze MGSIM results against ground truth.

    Parameters
    ----------
    results_path : str
        Path to combined MGSIM results NetCDF
    gt_path : str
        Path to ground truth CSV
    output_path : str, optional
        Path to save analysis NetCDF
    figure_path : str, optional
        Path to save figures
    """
    print("Loading data...")
    ds = xr.open_dataset(results_path)

    rows, cols = ds.dims['y'], ds.dims['x']
    ground_truth = load_ground_truth(gt_path, rows, cols)

    n_realizations = len(ds.realization)
    print(f"Loaded {n_realizations} realizations, grid: {rows} x {cols}")

    # Compute error for mean field
    print("\nComputing error for mean field...")
    mean_field = ds['mean'].values
    mean_metrics = compute_metrics(mean_field, ground_truth)

    print(f"  RMSE: {mean_metrics['rmse']:.4f}")
    print(f"  MAE:  {mean_metrics['mae']:.4f}")
    print(f"  R²:   {mean_metrics['r2']:.4f}")
    print(f"  Bias: {mean_metrics['bias']:.4f}")
    print(f"  SSIM: {mean_metrics['ssim']:.4f}")

    # Compute error for each realization
    print("\nComputing per-realization metrics...")
    realization_metrics = {
        'rmse': [],
        'mae': [],
        'r2': [],
        'bias': [],
        'ssim': [],
    }

    for i in range(n_realizations):
        realization = ds['simulated'].isel(realization=i).values
        metrics = compute_metrics(realization, ground_truth)
        for key in realization_metrics:
            realization_metrics[key].append(metrics[key])

        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{n_realizations} realizations")

    # Convert to arrays
    for key in realization_metrics:
        realization_metrics[key] = np.array(realization_metrics[key])

    # Compute error maps
    print("\nComputing error maps...")
    error_mean = mean_field - ground_truth
    abs_error_mean = np.abs(error_mean)

    # Create output dataset
    print("\nCreating output dataset...")
    ds_out = xr.Dataset(
        data_vars={
            'ground_truth': (['y', 'x'], ground_truth),
            'mean_prediction': (['y', 'x'], mean_field),
            'error': (['y', 'x'], error_mean),
            'abs_error': (['y', 'x'], abs_error_mean),
            'prediction_variance': (['y', 'x'], ds['variance'].values),
            'prediction_std': (['y', 'x'], ds['std'].values),
            # Per-realization metrics
            'realization_rmse': (['realization'], realization_metrics['rmse']),
            'realization_mae': (['realization'], realization_metrics['mae']),
            'realization_r2': (['realization'], realization_metrics['r2']),
            'realization_bias': (['realization'], realization_metrics['bias']),
            'realization_ssim': (['realization'], realization_metrics['ssim']),
        },
        coords={
            'y': ds.y,
            'x': ds.x,
            'realization': ds.realization,
        },
        attrs={
            'description': 'Error analysis of MGSIM results vs ground truth',
            'n_realizations': n_realizations,
            'mean_rmse': float(mean_metrics['rmse']),
            'mean_mae': float(mean_metrics['mae']),
            'mean_r2': float(mean_metrics['r2']),
            'mean_bias': float(mean_metrics['bias']),
            'mean_ssim': float(mean_metrics['ssim']),
        }
    )

    # Save output
    if output_path:
        print(f"\nSaving analysis to {output_path}...")
        ds_out.to_netcdf(output_path)

    # Generate figures
    if figure_path:
        print(f"\nGenerating figures...")
        generate_figures(ds_out, ground_truth, mean_field, figure_path, realization_metrics)

    # Print summary
    print("\n" + "=" * 70)
    print("ERROR ANALYSIS SUMMARY")
    print("=" * 70)
    print("\nMean Field Metrics:")
    print(f"  {'Metric':<12} {'Value':<15}")
    print(f"  {'-'*27}")
    for key, val in mean_metrics.items():
        print(f"  {key.upper():<12} {val:<15.4f}")

    print("\nPer-Realization Metrics (mean ± std across realizations):")
    print(f"  {'Metric':<12} {'Mean':<12} {'Std':<12} {'Min':<12} {'Max':<12}")
    print(f"  {'-'*60}")
    for key in ['rmse', 'mae', 'r2', 'bias', 'ssim']:
        arr = realization_metrics[key]
        print(f"  {key.upper():<12} {arr.mean():<12.4f} {arr.std():<12.4f} {arr.min():<12.4f} {arr.max():<12.4f}")

    return ds_out


def generate_figures(ds: xr.Dataset, ground_truth: np.ndarray, mean_field: np.ndarray,
                     output_dir: str, realization_metrics: dict):
    """Generate error analysis figures with geosoft colormap and hillshading."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Initialize colormap and light source
    geosoft_cmap = geosoft_cmap_k65()
    ls = LightSource(azdeg=315, altdeg=45)

    # Fixed color scale
    vmin, vmax = -1000, 1500
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Get extent from coordinates
    x_coords = ds.x.values
    y_coords = ds.y.values
    extent = [x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()]

    # Figure 1: Error maps with hillshading
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # Ground truth with hillshading
    ax = axes[0, 0]
    gt_rgb = ls.shade(ground_truth, cmap=geosoft_cmap, blend_mode='soft', vmin=vmin, vmax=vmax)
    ax.imshow(gt_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax.set_title('Ground Truth')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    sm = plt.cm.ScalarMappable(norm=norm, cmap=geosoft_cmap)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=0.6, label='Value')

    # Mean prediction with hillshading
    ax = axes[0, 1]
    mean_rgb = ls.shade(mean_field, cmap=geosoft_cmap, blend_mode='soft', vmin=vmin, vmax=vmax)
    ax.imshow(mean_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax.set_title(f'MGSIM Mean (R²={ds.attrs["mean_r2"]:.3f})')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(sm, ax=ax, shrink=0.6, label='Value')

    # Error (use diverging colormap, no hillshade)
    ax = axes[0, 2]
    error = ds['error'].values
    err_vmax = np.percentile(np.abs(error), 99)
    im = ax.imshow(error, cmap='RdBu_r', origin='lower', extent=extent,
                   vmin=-err_vmax, vmax=err_vmax, interpolation='nearest')
    ax.set_title(f'Error (Pred - True)\nRMSE={ds.attrs["mean_rmse"]:.3f}')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, shrink=0.6, label='Error')

    # Prediction variance with hillshading
    ax = axes[1, 0]
    variance_data = ds['prediction_variance'].values
    var_vmin, var_vmax = np.nanmin(variance_data), np.nanmax(variance_data)
    var_norm = Normalize(vmin=var_vmin, vmax=var_vmax)
    var_rgb = ls.shade(variance_data, cmap=geosoft_cmap, blend_mode='soft',
                       vmin=var_vmin, vmax=var_vmax)
    ax.imshow(var_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax.set_title('Prediction Variance')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    sm_var = plt.cm.ScalarMappable(norm=var_norm, cmap=geosoft_cmap)
    sm_var.set_array([])
    plt.colorbar(sm_var, ax=ax, shrink=0.6, label='Variance')

    # Absolute error with hillshading
    ax = axes[1, 1]
    abs_err = ds['abs_error'].values
    abs_vmin, abs_vmax = np.nanmin(abs_err), np.nanmax(abs_err)
    abs_norm = Normalize(vmin=abs_vmin, vmax=abs_vmax)
    abs_rgb = ls.shade(abs_err, cmap=plt.cm.Reds, blend_mode='soft',
                       vmin=abs_vmin, vmax=abs_vmax)
    ax.imshow(abs_rgb, origin='lower', extent=extent, interpolation='nearest')
    ax.set_title(f'Absolute Error\nMAE={ds.attrs["mean_mae"]:.3f}')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    sm_abs = plt.cm.ScalarMappable(norm=abs_norm, cmap=plt.cm.Reds)
    sm_abs.set_array([])
    plt.colorbar(sm_abs, ax=ax, shrink=0.6, label='|Error|')

    # Scatter: variance vs abs error
    ax = axes[1, 2]
    variance = ds['prediction_variance'].values.flatten()
    abs_error = ds['abs_error'].values.flatten()
    ax.scatter(variance, abs_error, alpha=0.1, s=1, c='#1f77b4')
    ax.set_xlabel('Prediction Variance')
    ax.set_ylabel('Absolute Error')
    ax.set_title('Variance vs Error')
    corr = np.corrcoef(variance, abs_error)[0, 1]
    ax.text(0.05, 0.95, f'Corr: {corr:.3f}', transform=ax.transAxes, va='top',
            fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_path / 'error_maps.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Figure 2: Metric distributions
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    metrics = ['rmse', 'mae', 'r2', 'bias', 'ssim']
    for i, metric in enumerate(metrics):
        ax = axes[i]
        data = realization_metrics[metric]
        ax.hist(data, bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(data.mean(), color='red', linestyle='--', label=f'Mean: {data.mean():.4f}')
        ax.set_xlabel(metric.upper())
        ax.set_ylabel('Count')
        ax.set_title(f'{metric.upper()} Distribution')
        ax.legend()

    axes[5].axis('off')

    plt.tight_layout()
    plt.savefig(output_path / 'metric_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved figures to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Compute error metrics against ground truth')
    parser.add_argument('--results', type=str, required=True,
                        help='Path to combined MGSIM results NetCDF')
    parser.add_argument('--ground-truth', type=str, required=True,
                        help='Path to ground truth CSV (x, y, val, clust)')
    parser.add_argument('--output', type=str, default='error_analysis.nc',
                        help='Output NetCDF file path')
    parser.add_argument('--figures', type=str, default='figures',
                        help='Directory to save figures')

    args = parser.parse_args()

    analyze_realizations(
        results_path=args.results,
        gt_path=args.ground_truth,
        output_path=args.output,
        figure_path=args.figures
    )


if __name__ == '__main__':
    main()
