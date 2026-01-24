import numpy as np
import pandas as pd
from typing import Optional


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
