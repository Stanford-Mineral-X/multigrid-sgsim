#!/usr/bin/env python
"""
Combine and Compare All 6 Ensembles

This script:
1. Combines partial NetCDF files for each ensemble
2. Computes error metrics vs ground truth for each
3. Generates comparison figures and tables

Usage:
    python combine_and_compare_ensembles.py --results-dir results --ground-truth gt_xyvc.csv
"""

import argparse
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from skimage.metrics import structural_similarity as ssim
from matplotlib.colors import LightSource


ENSEMBLE_NAMES = [
    'iso_subregions_dense',
    'iso_global_dense',
    'aniso_subregions_dense',
    'aniso_global_dense',
    'iso_subregions_medium',
    'iso_subregions_sparse',
]

ENSEMBLE_LABELS = {
    'iso_subregions_dense': 'Iso + Subregions (Dense)',
    'iso_global_dense': 'Iso + Global (Dense)',
    'aniso_subregions_dense': 'Aniso + Subregions (Dense)',
    'aniso_global_dense': 'Aniso + Global (Dense)',
    'iso_subregions_medium': 'Iso + Subregions (Medium)',
    'iso_subregions_sparse': 'Iso + Subregions (Sparse)',
}


def combine_ensemble(input_dir: Path, ensemble_name: str) -> xr.Dataset:
    """Combine partial files for one ensemble."""
    ensemble_dir = input_dir / ensemble_name
    files = sorted(ensemble_dir.glob('realizations_*.nc'))

    if not files:
        print(f"  No files found for {ensemble_name}")
        return None

    print(f"  Found {len(files)} files")

    datasets = []
    for f in files:
        ds = xr.open_dataset(f)
        datasets.append(ds)

    combined = xr.concat(datasets, dim='realization')
    combined = combined.sortby('realization')

    # Compute statistics
    simulated = combined['simulated']
    combined['mean'] = simulated.mean(dim='realization')
    combined['variance'] = simulated.var(dim='realization')
    combined['std'] = simulated.std(dim='realization')
    combined['median'] = simulated.median(dim='realization')

    combined.attrs['ensemble_name'] = ensemble_name
    combined.attrs['n_realizations'] = len(combined.realization)

    return combined


def compute_metrics(pred: np.ndarray, true: np.ndarray) -> dict:
    """Compute error metrics."""
    pred_flat = pred.flatten()
    true_flat = true.flatten()

    mask = ~(np.isnan(pred_flat) | np.isnan(true_flat))
    pred_valid = pred_flat[mask]
    true_valid = true_flat[mask]

    rmse = np.sqrt(mean_squared_error(true_valid, pred_valid))
    mae = mean_absolute_error(true_valid, pred_valid)
    r2 = r2_score(true_valid, pred_valid)
    bias = np.mean(pred_valid - true_valid)

    pred_norm = (pred - np.nanmin(pred)) / (np.nanmax(pred) - np.nanmin(pred) + 1e-10)
    true_norm = (true - np.nanmin(true)) / (np.nanmax(true) - np.nanmin(true) + 1e-10)
    ssim_val = ssim(true_norm, pred_norm, data_range=1.0)

    return {'rmse': rmse, 'mae': mae, 'r2': r2, 'bias': bias, 'ssim': ssim_val}


def analyze_ensemble(ds: xr.Dataset, ground_truth: np.ndarray) -> dict:
    """Analyze one ensemble against ground truth."""
    mean_field = ds['mean'].values

    # Mean field metrics
    mean_metrics = compute_metrics(mean_field, ground_truth)

    # Per-realization metrics
    n_real = len(ds.realization)
    real_metrics = {k: [] for k in ['rmse', 'mae', 'r2', 'bias', 'ssim']}

    for i in range(n_real):
        realization = ds['simulated'].isel(realization=i).values
        metrics = compute_metrics(realization, ground_truth)
        for k, v in metrics.items():
            real_metrics[k].append(v)

    for k in real_metrics:
        real_metrics[k] = np.array(real_metrics[k])

    return {
        'mean_metrics': mean_metrics,
        'realization_metrics': real_metrics,
        'error_field': mean_field - ground_truth,
        'variance_field': ds['variance'].values,
    }


def generate_comparison_figures(results: dict, ground_truth: np.ndarray, output_dir: Path):
    """Generate comparison figures."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Figure 1: Error maps for all ensembles (2 rows: error, variance)
    n_ensembles = len([n for n in ENSEMBLE_NAMES if n in results])
    fig, axes = plt.subplots(2, n_ensembles, figsize=(5 * n_ensembles, 10))

    # Handle case where only 1 ensemble
    if n_ensembles == 1:
        axes = axes.reshape(2, 1)

    # Find common error scale
    all_errors = [results[name]['error_field'] for name in ENSEMBLE_NAMES if name in results]
    vmax = np.percentile(np.abs(np.concatenate([e.flatten() for e in all_errors])), 99)

    col = 0
    for name in ENSEMBLE_NAMES:
        if name not in results:
            continue
        res = results[name]

        # Error map
        ax = axes[0, col]
        im = ax.imshow(res['error_field'], cmap='RdBu_r', origin='lower', vmin=-vmax, vmax=vmax)
        ax.set_title(f"{ENSEMBLE_LABELS[name]}\nRMSE={res['mean_metrics']['rmse']:.3f}")
        plt.colorbar(im, ax=ax, label='Error')

        # Variance map
        ax = axes[1, col]
        im = ax.imshow(res['variance_field'], cmap='YlOrRd', origin='lower')
        ax.set_title(f"Prediction Variance\nR²={res['mean_metrics']['r2']:.3f}")
        plt.colorbar(im, ax=ax, label='Variance')

        col += 1

    plt.suptitle(f'Error Maps: {n_ensembles} Ensemble Comparison', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'ensemble_error_maps.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Figure 2: Metric distributions as box plots
    fig, axes = plt.subplots(1, 5, figsize=(20, 5))

    metrics = ['rmse', 'mae', 'r2', 'bias', 'ssim']

    for i, metric in enumerate(metrics):
        ax = axes[i]
        data = []
        labels = []
        for name in ENSEMBLE_NAMES:
            if name in results:
                data.append(results[name]['realization_metrics'][metric])
                labels.append(ENSEMBLE_LABELS[name].replace(' + ', '\n'))

        bp = ax.boxplot(data, labels=labels, patch_artist=True)
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']  # 6 colors
        for patch, color in zip(bp['boxes'], colors[:len(data)]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel(metric.upper())
        ax.set_title(f'{metric.upper()} Distribution')
        ax.tick_params(axis='x', rotation=45)

    plt.suptitle('Per-Realization Metric Distributions', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_dir / 'ensemble_metric_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Figure 3: Bar chart comparing mean metrics
    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(metrics))
    width = 0.12  # Narrower for 6 ensembles
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']  # 6 colors

    for i, name in enumerate(ENSEMBLE_NAMES):
        if name not in results:
            continue
        values = [results[name]['mean_metrics'][m] for m in metrics]
        # Normalize for visualization (different scales)
        ax.bar(x + i * width, values, width, label=ENSEMBLE_LABELS[name], color=colors[i], alpha=0.8)

    ax.set_ylabel('Metric Value')
    ax.set_title('Mean Field Metrics Comparison')
    ax.set_xticks(x + width * 2.5)  # Center for 6 bars
    ax.set_xticklabels([m.upper() for m in metrics])
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'ensemble_metric_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Figures saved to {output_dir}")


def generate_hillshade_figures(combined_datasets: dict, output_dir: Path, cmap='viridis'):
    """Generate hillshade figures showing ensemble mean + 3 random realizations.

    Creates one stacked figure (1 column, 4 rows) per ensemble.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    ls = LightSource(azdeg=315, altdeg=45)

    for name, ds in combined_datasets.items():
        n_real = len(ds.realization)
        if n_real < 3:
            print(f"  Skipping {name}: need at least 3 realizations, have {n_real}")
            continue

        # Get ensemble mean
        mean_field = ds['mean'].values

        # Select 3 random realizations
        np.random.seed(42)  # For reproducibility
        random_indices = np.random.choice(n_real, size=3, replace=False)
        random_indices.sort()

        # Compute common color scale across all 4 panels
        all_fields = [mean_field] + [ds['simulated'].isel(realization=i).values for i in random_indices]
        vmin = np.nanpercentile(np.concatenate([f.flatten() for f in all_fields]), 2)
        vmax = np.nanpercentile(np.concatenate([f.flatten() for f in all_fields]), 98)

        # Create figure: 1 column, 4 rows
        fig, axes = plt.subplots(4, 1, figsize=(10, 24))

        # Get colormap
        if isinstance(cmap, str):
            cmap_obj = plt.get_cmap(cmap)
        else:
            cmap_obj = cmap

        # Panel 0: Ensemble mean
        ax = axes[0]
        mean_rgb = ls.shade(mean_field, cmap=cmap_obj, blend_mode='soft', vmin=vmin, vmax=vmax)
        ax.imshow(mean_rgb, origin='lower', interpolation='nearest')
        ax.set_title(f'{ENSEMBLE_LABELS[name]}\nEnsemble Mean (n={n_real})', fontsize=12)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        # Panels 1-3: Random realizations
        for i, real_idx in enumerate(random_indices):
            ax = axes[i + 1]
            realization = ds['simulated'].isel(realization=real_idx).values
            real_rgb = ls.shade(realization, cmap=cmap_obj, blend_mode='soft', vmin=vmin, vmax=vmax)
            ax.imshow(real_rgb, origin='lower', interpolation='nearest')
            ax.set_title(f'Realization {int(real_idx)}', fontsize=12)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Value')

        plt.tight_layout()
        plt.savefig(output_dir / f'{name}_hillshade_realizations.png', dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  Saved hillshade figure for {name}")


def generate_latex_table(results: dict) -> str:
    """Generate LaTeX table for paper."""
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Comparison of MGSIM ensemble configurations against ground truth}",
        r"\label{tab:ensemble_comparison}",
        r"\begin{tabular}{lcccccc}",
        r"\hline",
        r"Configuration & RMSE & MAE & R$^2$ & Bias & SSIM \\",
        r"\hline",
    ]

    for name in ENSEMBLE_NAMES:
        if name not in results:
            continue
        m = results[name]['mean_metrics']
        label = ENSEMBLE_LABELS[name]
        lines.append(
            f"{label} & {m['rmse']:.4f} & {m['mae']:.4f} & {m['r2']:.4f} & {m['bias']:.4f} & {m['ssim']:.4f} \\\\"
        )

    lines.extend([
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
    ])

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description='Combine and compare all ensembles')
    parser.add_argument('--results-dir', type=str, default='results',
                        help='Directory containing ensemble subdirectories')
    parser.add_argument('--ground-truth', type=str, required=True,
                        help='Path to ground truth CSV')
    parser.add_argument('--output-dir', type=str, default='analysis',
                        help='Output directory for figures and tables')
    parser.add_argument('--cmap', type=str, default='viridis',
                        help='Colormap for hillshade figures (default: viridis)')

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load ground truth
    print("Loading ground truth...")
    df_gt = pd.read_csv(args.ground_truth)
    # Determine grid shape from data
    x_unique = df_gt['x'].nunique()
    y_unique = df_gt['y'].nunique()
    ground_truth = df_gt['val'].values.reshape(y_unique, x_unique)
    print(f"  Grid shape: {ground_truth.shape}")

    # Process each ensemble
    results = {}
    combined_datasets = {}

    for name in ENSEMBLE_NAMES:
        print(f"\nProcessing ensemble: {name}")

        # Combine partial files
        ds = combine_ensemble(results_dir, name)
        if ds is None:
            continue

        combined_datasets[name] = ds
        print(f"  Combined {len(ds.realization)} realizations")

        # Save combined file
        combined_path = output_dir / f'{name}_combined.nc'
        ds.to_netcdf(combined_path)
        print(f"  Saved to {combined_path}")

        # Analyze
        print("  Computing metrics...")
        results[name] = analyze_ensemble(ds, ground_truth)

    # Generate comparison figures
    print("\nGenerating comparison figures...")
    generate_comparison_figures(results, ground_truth, output_dir / 'figures')

    # Generate hillshade figures (mean + 3 random realizations per ensemble)
    print("\nGenerating hillshade realization figures...")
    generate_hillshade_figures(combined_datasets, output_dir / 'figures', cmap=args.cmap)

    # Generate LaTeX table
    latex_table = generate_latex_table(results)
    with open(output_dir / 'comparison_table.tex', 'w') as f:
        f.write(latex_table)
    print(f"\nLaTeX table saved to {output_dir / 'comparison_table.tex'}")

    # Print summary
    print("\n" + "=" * 80)
    print("ENSEMBLE COMPARISON SUMMARY")
    print("=" * 80)
    print(f"\n{'Configuration':<25} {'RMSE':<10} {'MAE':<10} {'R²':<10} {'SSIM':<10}")
    print("-" * 65)

    for name in ENSEMBLE_NAMES:
        if name not in results:
            continue
        m = results[name]['mean_metrics']
        print(f"{ENSEMBLE_LABELS[name]:<25} {m['rmse']:<10.4f} {m['mae']:<10.4f} {m['r2']:<10.4f} {m['ssim']:<10.4f}")

    # Highlight best performer
    print("\n" + "-" * 65)
    best_rmse = min(ENSEMBLE_NAMES, key=lambda n: results[n]['mean_metrics']['rmse'] if n in results else float('inf'))
    best_r2 = max(ENSEMBLE_NAMES, key=lambda n: results[n]['mean_metrics']['r2'] if n in results else float('-inf'))
    print(f"Best RMSE: {ENSEMBLE_LABELS[best_rmse]}")
    print(f"Best R²:   {ENSEMBLE_LABELS[best_r2]}")

    # Save results as pickle for further analysis
    import pickle
    with open(output_dir / 'ensemble_results.pkl', 'wb') as f:
        pickle.dump(results, f)


if __name__ == '__main__':
    main()
