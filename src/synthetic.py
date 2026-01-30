import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple


def extract_flightlines(
    df: pd.DataFrame,
    angle: float = 0.0,
    spacing: int = 10,
    x_col: str = 'x',
    y_col: str = 'y',
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Extract synthetic flightline data from a gridded dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Input gridded data with x, y coordinates and values
    angle : float
        Angle of flightlines relative to x-axis in degrees.
        0 = horizontal (parallel to x), 90 = vertical (parallel to y)
    spacing : int
        Number of grid cells between flightlines along the x-axis
    x_col, y_col : str
        Column names for x and y coordinates
    save_path : str, optional
        If provided, save the extracted flightlines to this path as CSV

    Returns
    -------
    pd.DataFrame
        Subset of input data representing synthetic flightlines
    """
    df = df.copy()

    # Get grid parameters
    x_vals = np.sort(df[x_col].unique())
    y_vals = np.sort(df[y_col].unique())
    dx = x_vals[1] - x_vals[0] if len(x_vals) > 1 else 1.0

    # Convert angle to radians
    theta = np.radians(angle)

    # Spacing measured along x-axis (in coordinate units)
    spacing_dist = spacing * dx

    # For a line at angle theta passing through x-intercept x0:
    #   y = tan(theta) * (x - x0)
    # Rearranged: x0 = x - y / tan(theta)
    #
    # We want lines spaced by `spacing_dist` along the x-axis,
    # so x-intercepts are at 0, spacing_dist, 2*spacing_dist, etc.

    if np.abs(np.cos(theta)) < 1e-10:
        # Near-vertical lines: use x directly
        df['_x_intercept'] = df[x_col]
    else:
        # Compute where this point's flightline crosses y=0
        df['_x_intercept'] = df[x_col] - df[y_col] / np.tan(theta)

    # Assign each point to nearest flightline
    df['_line_idx'] = np.round(df['_x_intercept'] / spacing_dist)

    # Keep points that fall on flightlines (within tolerance)
    tolerance = dx * 0.5
    df['_dist_to_line'] = np.abs(df['_x_intercept'] - df['_line_idx'] * spacing_dist)

    flightlines = df[df['_dist_to_line'] < tolerance].copy()

    # Clean up temporary columns
    flightlines = flightlines.drop(columns=['_x_intercept', '_line_idx', '_dist_to_line'])

    # Save if requested
    if save_path is not None:
        flightlines.to_csv(save_path, index=False)
        print(f"Saved {len(flightlines)} flightline points to {save_path}")

    return flightlines


def extract_nested_flightlines(
    df: pd.DataFrame,
    spacings: List[int],
    angle: float = 0.0,
    x_col: str = 'x',
    y_col: str = 'y',
    save_dir: Optional[str] = None,
    names: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Extract multiple nested flightline datasets where sparser datasets
    are guaranteed to be subsets of denser datasets.

    Parameters
    ----------
    df : pd.DataFrame
        Input gridded data with x, y coordinates and values
    spacings : List[int]
        List of spacings in ascending order (densest first).
        Each spacing should be a multiple of the previous one.
        Example: [3, 6, 12] means lines every 3, 6, and 12 cells.
    angle : float
        Angle of flightlines relative to x-axis in degrees.
        0 = horizontal (parallel to x), 90 = vertical (parallel to y)
    x_col, y_col : str
        Column names for x and y coordinates
    save_dir : str, optional
        If provided, save each dataset to this directory
    names : List[str], optional
        Names for each spacing level (for keys and filenames).
        Defaults to ['spacing_3', 'spacing_6', etc.]

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping names to flightline DataFrames

    Example
    -------
    >>> spacings = [3, 6]  # dense=3, sparse=6
    >>> datasets = extract_nested_flightlines(df_gt, spacings, angle=0)
    >>> df_dense = datasets['spacing_3']   # 33 flight lines
    >>> df_sparse = datasets['spacing_6']  # 17 flight lines (subset of dense)
    """
    # Validate spacings are in ascending order and nested
    spacings = sorted(spacings)
    base_spacing = spacings[0]
    for s in spacings[1:]:
        if s % base_spacing != 0:
            raise ValueError(
                f"Spacing {s} is not a multiple of base spacing {base_spacing}. "
                f"All spacings must be multiples of the smallest spacing to ensure nesting."
            )

    # Generate names if not provided
    if names is None:
        names = [f"spacing_{s}" for s in spacings]
    elif len(names) != len(spacings):
        raise ValueError(f"Length of names ({len(names)}) must match length of spacings ({len(spacings)})")

    # Extract densest dataset first
    df_dense = extract_flightlines(df, angle=angle, spacing=base_spacing, x_col=x_col, y_col=y_col)

    # Get the line indices from the dense dataset
    x_vals = np.sort(df[x_col].unique())
    dx = x_vals[1] - x_vals[0] if len(x_vals) > 1 else 1.0
    theta = np.radians(angle)

    if np.abs(np.cos(theta)) < 1e-10:
        df_dense['_line_coord'] = df_dense[x_col]
    else:
        df_dense['_line_coord'] = df_dense[x_col] - df_dense[y_col] / np.tan(theta)

    # Round to nearest line based on base spacing
    base_spacing_dist = base_spacing * dx
    df_dense['_line_idx'] = np.round(df_dense['_line_coord'] / base_spacing_dist).astype(int)

    results = {}

    for spacing, name in zip(spacings, names):
        # Calculate which line indices to keep (multiples of spacing/base_spacing)
        step = spacing // base_spacing

        # Keep lines where line_idx is a multiple of step
        mask = (df_dense['_line_idx'] % step) == 0
        df_subset = df_dense[mask].copy()

        # Clean up temporary columns
        df_subset = df_subset.drop(columns=['_line_coord', '_line_idx'], errors='ignore')

        # Count unique lines
        if np.abs(np.cos(theta)) < 1e-10:
            n_lines = df_subset[x_col].nunique()
        else:
            line_coords = df_subset[x_col] - df_subset[y_col] / np.tan(theta)
            n_lines = len(np.unique(np.round(line_coords / (spacing * dx))))

        print(f"  {name}: {len(df_subset)} points on {n_lines} flight lines (spacing={spacing})")

        # Save if requested
        if save_dir is not None:
            from pathlib import Path
            save_path = Path(save_dir) / f"fl_{name}.csv"
            df_subset.to_csv(save_path, index=False)
            print(f"    Saved to {save_path}")

        results[name] = df_subset

    # Verify nesting property
    print("\nVerifying nesting property:")
    for i in range(len(spacings) - 1):
        dense_name = names[i]
        sparse_name = names[i + 1]
        dense_set = set(zip(results[dense_name][x_col], results[dense_name][y_col]))
        sparse_set = set(zip(results[sparse_name][x_col], results[sparse_name][y_col]))

        if sparse_set.issubset(dense_set):
            print(f"  ✓ {sparse_name} is subset of {dense_name}")
        else:
            missing = sparse_set - dense_set
            print(f"  ✗ WARNING: {sparse_name} has {len(missing)} points not in {dense_name}")

    return results


def create_flightline_experiment(
    df_ground_truth: pd.DataFrame,
    dense_spacing: int = 3,
    sparse_spacing: int = 6,
    angle: float = 0.0,
    x_col: str = 'x',
    y_col: str = 'y',
    save_dir: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create dense and sparse flightline datasets for experimental comparison.

    This is a convenience function for the common case of comparing
    two flight line densities.

    Parameters
    ----------
    df_ground_truth : pd.DataFrame
        Ground truth gridded data
    dense_spacing : int
        Spacing for dense flight lines (in grid cells)
    sparse_spacing : int
        Spacing for sparse flight lines (must be multiple of dense_spacing)
    angle : float
        Flight line angle in degrees (0 = E-W, 90 = N-S)
    x_col, y_col : str
        Column names for coordinates in the ground truth DataFrame
    save_dir : str, optional
        Directory to save CSV files

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (df_dense, df_sparse) - The two flight line datasets

    Example
    -------
    >>> df_gt = pd.read_csv('gt_xyvc.csv')
    >>> df_dense, df_sparse = create_flightline_experiment(
    ...     df_gt, dense_spacing=3, sparse_spacing=6, save_dir='data/'
    ... )
    >>> print(f"Dense: {len(df_dense)} points, Sparse: {len(df_sparse)} points")
    """
    if sparse_spacing % dense_spacing != 0:
        raise ValueError(
            f"sparse_spacing ({sparse_spacing}) must be a multiple of "
            f"dense_spacing ({dense_spacing}) to ensure nesting."
        )

    print(f"Creating flight line experiment:")
    print(f"  Dense spacing: {dense_spacing} cells")
    print(f"  Sparse spacing: {sparse_spacing} cells")
    print(f"  Angle: {angle}°")

    datasets = extract_nested_flightlines(
        df_ground_truth,
        spacings=[dense_spacing, sparse_spacing],
        angle=angle,
        x_col=x_col,
        y_col=y_col,
        save_dir=save_dir,
        names=['dense', 'sparse'],
    )

    return datasets['dense'], datasets['sparse']
