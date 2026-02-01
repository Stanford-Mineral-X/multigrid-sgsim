#!/usr/bin/env python
"""
Generate Synthetic Flight Line Datasets for Ensemble Experiments

This script generates nested flight line datasets from ground truth data:
- Dense: flight lines every 4 grid cells (gap of 3)
- Sparse: flight lines every 8 grid cells (gap of 7, subset of dense)

The sparse dataset is guaranteed to be a subset of the dense dataset,
ensuring a fair comparison of flight line spacing effects.

Usage:
    python generate_flightline_datasets.py

Output:
    data/fl_dense.csv   - Dense flight lines (spacing=4, gap=3)
    data/fl_sparse.csv  - Sparse flight lines (spacing=8, gap=7)
"""

import sys
from pathlib import Path

# Add src directory to path
script_dir = Path(__file__).parent
src_dir = script_dir.parent / 'src'
sys.path.insert(0, str(src_dir))

import numpy as np
import pandas as pd
from synthetic import extract_nested_flightlines


# =============================================================================
# USER PARAMETERS - MODIFY THESE FOR YOUR DATASET
# =============================================================================

# Input ground truth file
GT_PATH = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data/gt_xyvc.csv')

# Output directory for flight line CSVs
OUTPUT_DIR = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/multigrid-sgsim/demos/data')

# Flight line parameters
DENSE_SPACING = 4   # Grid cells between dense flight lines (gap of 3)
SPARSE_SPACING = 8  # Grid cells between sparse flight lines (gap of 7, must be multiple of dense)
ANGLE = -45.0         # Flight line angle: 0 = E-W, 90 = N-S

# Column names in ground truth file
X_COL = 'x'
Y_COL = 'y'
VALUE_COL = 'val'
CLUSTER_COL = 'clust'


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("GENERATE SYNTHETIC FLIGHT LINE DATASETS")
    print("=" * 70)

    # Load ground truth
    print(f"\nLoading ground truth from: {GT_PATH}")
    df_gt = pd.read_csv(GT_PATH)

    # Standardize column names
    df_gt.columns = [X_COL, Y_COL, VALUE_COL, CLUSTER_COL]

    # Get grid info
    x_unique = np.sort(df_gt[X_COL].unique())
    y_unique = np.sort(df_gt[Y_COL].unique())
    dx = x_unique[1] - x_unique[0] if len(x_unique) > 1 else 1.0
    dy = y_unique[1] - y_unique[0] if len(y_unique) > 1 else 1.0

    print(f"  Grid size: {len(x_unique)} x {len(y_unique)}")
    print(f"  Grid spacing: dx={dx}, dy={dy}")
    print(f"  X range: [{x_unique.min()}, {x_unique.max()}]")
    print(f"  Y range: [{y_unique.min()}, {y_unique.max()}]")
    print(f"  Total points: {len(df_gt)}")
    print(f"  Clusters: {df_gt[CLUSTER_COL].nunique()}")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Generate flight line datasets
    print(f"\nGenerating flight line datasets:")
    print(f"  Dense spacing: {DENSE_SPACING} cells ({DENSE_SPACING * dx} units)")
    print(f"  Sparse spacing: {SPARSE_SPACING} cells ({SPARSE_SPACING * dx} units)")
    print(f"  Angle: {ANGLE}° (0°=E-W, 90°=N-S)")

    # Use extract_nested_flightlines directly for custom file naming
    datasets = extract_nested_flightlines(
        df_gt,
        spacings=[DENSE_SPACING, SPARSE_SPACING],
        angle=ANGLE,
        x_col=X_COL,
        y_col=Y_COL,
        save_dir=str(OUTPUT_DIR),
        names=['xyvc_dense', 'xyvc_sparse'],
    )
    df_dense = datasets['xyvc_dense']
    df_sparse = datasets['xyvc_sparse']

    # Compute statistics
    print("\n" + "=" * 70)
    print("DATASET SUMMARY")
    print("=" * 70)

    for name, df in [("Dense", df_dense), ("Sparse", df_sparse)]:
        n_points = len(df)
        n_lines = df[Y_COL].nunique() if ANGLE == 0 else df[X_COL].nunique()
        coverage = n_points / len(df_gt) * 100

        print(f"\n{name} Dataset:")
        print(f"  Points: {n_points}")
        print(f"  Flight lines: {n_lines}")
        print(f"  Coverage: {coverage:.1f}% of ground truth")
        print(f"  Value range: [{df[VALUE_COL].min():.2f}, {df[VALUE_COL].max():.2f}]")
        print(f"  Value mean: {df[VALUE_COL].mean():.2f}")
        print(f"  Value std: {df[VALUE_COL].std():.2f}")

        # Points per cluster
        print(f"  Points per cluster:")
        for c in sorted(df[CLUSTER_COL].unique()):
            n = len(df[df[CLUSTER_COL] == c])
            print(f"    Cluster {int(c)}: {n} points")

    # Verify nesting
    dense_coords = set(zip(df_dense[X_COL], df_dense[Y_COL]))
    sparse_coords = set(zip(df_sparse[X_COL], df_sparse[Y_COL]))

    print("\n" + "=" * 70)
    print("NESTING VERIFICATION")
    print("=" * 70)
    if sparse_coords.issubset(dense_coords):
        print("✓ PASSED: Sparse dataset is a proper subset of dense dataset")
        print(f"  Dense has {len(dense_coords)} unique locations")
        print(f"  Sparse has {len(sparse_coords)} unique locations")
        print(f"  Sparse/Dense ratio: {len(sparse_coords)/len(dense_coords)*100:.1f}%")
    else:
        missing = sparse_coords - dense_coords
        print(f"✗ FAILED: Sparse has {len(missing)} points not in dense!")

    # Print output file locations
    print("\n" + "=" * 70)
    print("OUTPUT FILES")
    print("=" * 70)
    print(f"  Dense:  {OUTPUT_DIR / 'fl_xyvc_dense.csv'}")
    print(f"  Sparse: {OUTPUT_DIR / 'fl_xyvc_sparse.csv'}")

    print("\nDone!")


if __name__ == '__main__':
    main()
