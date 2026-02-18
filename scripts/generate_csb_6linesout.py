#!/usr/bin/env python
"""
Generate CSB holdout dataset by withholding 6 N-S flight lines.

Produces:
  - csb_xyvmmc_6linesout.csv : same format as original, fl_mask=0 for withheld points
  - csb_holdout_6lines.csv   : just the withheld observations (for validation)
  - csb_holdout_6lines.json  : metadata (x values, tolerance, counts)

Mirrors the holdout logic from sandbox/csb_work.ipynb.

Usage:
    python generate_csb_6linesout.py
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path


# =============================================================================
# USER PARAMETERS
# =============================================================================

DATA_DIR = Path('/Users/jrines/stanford_gp/research/mx/computers_geosciences/data')
CSB_PATH = DATA_DIR / 'csb_xyvmmc.csv'

# Flight lines to hold out (x-coordinates of N-S lines)
HOLDOUT_X_VALUES = [548400, 551000, 553650, 555050, 558610, 560950]
TOL = 55  # tolerance in coordinate units


def main():
    print("=" * 60)
    print("CSB HOLDOUT DATASET GENERATOR — 6 LINES")
    print("=" * 60)

    # Load original data
    df = pd.read_csv(CSB_PATH)
    n_total = len(df)
    n_fl_orig = (df['fl_mask'] == 1).sum()
    print(f"Original: {n_total} grid points, {n_fl_orig} flight line observations")

    # Identify flight line points to hold out
    fl_mask = df['fl_mask'] == 1
    holdout_mask = pd.Series(False, index=df.index)

    for x_val in HOLDOUT_X_VALUES:
        holdout_mask |= fl_mask & ((df['x'] - x_val).abs() <= TOL)

    n_holdout = holdout_mask.sum()
    print(f"\nHolding out {n_holdout} points at x = {HOLDOUT_X_VALUES} (tol={TOL})")

    # Save holdout observations (ground truth for validation)
    df_holdout = df.loc[holdout_mask, ['x', 'y', 'val', 'cluster']].copy()
    holdout_path = DATA_DIR / 'csb_holdout_6lines.csv'
    df_holdout.to_csv(holdout_path, index=False)
    print(f"Saved holdout data: {holdout_path} ({len(df_holdout)} points)")

    # Create reduced dataset: set fl_mask=0 for withheld points
    df_reduced = df.copy()
    df_reduced.loc[holdout_mask, 'fl_mask'] = 0
    n_fl_reduced = (df_reduced['fl_mask'] == 1).sum()

    reduced_path = DATA_DIR / 'csb_xyvmmc_6linesout.csv'
    df_reduced.to_csv(reduced_path, index=False)
    print(f"Saved reduced dataset: {reduced_path} ({n_fl_reduced} flight line obs)")

    # Save metadata
    per_line = {}
    for x_val in HOLDOUT_X_VALUES:
        line_mask = fl_mask & ((df['x'] - x_val).abs() <= TOL)
        per_line[str(x_val)] = {
            'n_points': int(line_mask.sum()),
            'y_min': float(df.loc[line_mask, 'y'].min()) if line_mask.any() else None,
            'y_max': float(df.loc[line_mask, 'y'].max()) if line_mask.any() else None,
        }

    metadata = {
        'holdout_x_values': HOLDOUT_X_VALUES,
        'tolerance': TOL,
        'n_holdout_total': int(n_holdout),
        'n_fl_original': int(n_fl_orig),
        'n_fl_reduced': int(n_fl_reduced),
        'per_line': per_line,
        'holdout_file': 'csb_holdout_6lines.csv',
        'reduced_file': 'csb_xyvmmc_6linesout.csv',
    }

    meta_path = DATA_DIR / 'csb_holdout_6lines.json'
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"Saved metadata: {meta_path}")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Original flight line obs:  {n_fl_orig}")
    print(f"  Held out:                  {n_holdout}")
    print(f"  Retained flight line obs:  {n_fl_reduced}")
    for x_val, info in per_line.items():
        if info['n_points'] > 0:
            print(f"    x={x_val}: {info['n_points']} points "
                  f"(y: {info['y_min']:.0f} to {info['y_max']:.0f})")
        else:
            print(f"    x={x_val}: 0 points (no flight line at this x)")


if __name__ == '__main__':
    main()
