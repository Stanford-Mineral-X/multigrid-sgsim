"""
Validation module for multigrid-sgsim package.

Provides tools for ensemble falsification testing using Robust Mahalanobis
Distance (RMD) and multi-method comparison visualization.

Functions
---------
robust_mahalanobis_distance
    Compute RMD for ensemble falsification testing.
assign_line_ids
    Assign line IDs to points using KMeans clustering.
extract_ensemble_along_line
    Extract ensemble values along a validation line.
create_validation_animation
    Create multi-method comparison animation with RMD.
plot_validation_frame
    Plot single validation frame for static use.

Classes
-------
ValidationAnimationConfig
    Configuration dataclass for animation parameters.

Notes
-----
RMD implementation based on David Yin's RobustMD_flsification.py (April 2019).
Contact: yinzhen@stanford.edu
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation, FFMpegWriter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Tuple
from scipy import stats
from sklearn.covariance import MinCovDet
from sklearn.cluster import KMeans

from utils import geosoft_cmap_k65


# =============================================================================
# RMD Falsification
# =============================================================================

def robust_mahalanobis_distance(
    ensemble: np.ndarray | xr.DataArray,
    observation: np.ndarray,
    quantile: float = 95.0,
    random_state: int | None = 0,
) -> tuple[float, float, np.ndarray]:
    """
    Compute Robust Mahalanobis Distance for ensemble falsification testing.

    Uses sklearn's MinCovDet for robust covariance estimation. Tests whether
    an observation is statistically consistent with an ensemble of realizations.

    Parameters
    ----------
    ensemble : np.ndarray or xr.DataArray
        Ensemble realizations with shape (n_realizations, n_points).
        If xr.DataArray, values are extracted automatically.
    observation : np.ndarray
        Observed values with shape (n_points,) or (1, n_points).
    quantile : float, default=95.0
        Percentile threshold for falsification (e.g., 95 or 97.5).
    random_state : int or None, default=0
        Random seed for MinCovDet reproducibility.

    Returns
    -------
    rmd_observation : float
        RMD of the observation relative to the ensemble.
    rmd_threshold : float
        The quantile threshold value from the ensemble RMD distribution.
    rmd_realizations : np.ndarray
        RMD values for each realization, shape (n_realizations,).

    Notes
    -----
    Falsification criterion: If rmd_observation > rmd_threshold, the
    observation is inconsistent with the ensemble at the given quantile level.

    Implementation based on David Yin's RobustMD_flsification.py (April 2019).

    Examples
    --------
    >>> ensemble = np.random.randn(1000, 50)  # 1000 realizations, 50 points
    >>> observation = np.random.randn(50)
    >>> rmd_obs, threshold, rmd_all = robust_mahalanobis_distance(ensemble, observation)
    >>> is_falsified = rmd_obs > threshold
    """
    # Handle xarray input
    if isinstance(ensemble, xr.DataArray):
        ensemble = ensemble.values

    # Ensure 2D array
    ensemble = np.atleast_2d(ensemble)
    n_realizations, n_points = ensemble.shape

    # Reshape observation to (1, n_points)
    observation = np.atleast_2d(observation)
    if observation.shape[0] != 1:
        observation = observation.T
    if observation.shape[1] != n_points:
        raise ValueError(
            f"Observation has {observation.shape[1]} points but ensemble has {n_points}"
        )

    # Fit robust covariance estimator
    mcd = MinCovDet(random_state=random_state).fit(ensemble)

    # Compute RMD for observation
    obs_centered = observation - mcd.location_
    cov_inv = np.linalg.inv(mcd.covariance_)
    rmd_observation = float(np.sqrt(obs_centered @ cov_inv @ obs_centered.T))

    # Compute RMD for each realization
    rmd_realizations = np.zeros(n_realizations)
    for i in range(n_realizations):
        sample_centered = ensemble[i:i+1, :] - mcd.location_
        rmd_realizations[i] = np.sqrt(sample_centered @ cov_inv @ sample_centered.T)

    # Compute quantile threshold
    rmd_threshold = float(stats.scoreatpercentile(rmd_realizations, quantile))

    return rmd_observation, rmd_threshold, rmd_realizations


# =============================================================================
# Line Extraction Utilities
# =============================================================================

def assign_line_ids(
    df: pd.DataFrame,
    angle_deg: float,
    spacing: float,
    x_col: str = 'x',
    y_col: str = 'y',
) -> pd.DataFrame:
    """
    Assign line IDs to points using KMeans clustering on perpendicular axis.

    Projects points onto the axis perpendicular to the flight line direction,
    then clusters them to identify distinct lines.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with x, y coordinates.
    angle_deg : float
        Flight line angle in degrees (0 = E-W, 90 = N-S, -45 = diagonal SW-NE).
    spacing : float
        Approximate spacing between lines (used to estimate n_clusters).
    x_col, y_col : str
        Column names for coordinates.

    Returns
    -------
    pd.DataFrame
        Copy of input with added 'line_id' column (spatially ordered).

    Notes
    -----
    The number of clusters is estimated as:
        n_clusters = round(perpendicular_extent / spacing) + 1

    Line IDs are remapped so they increase spatially along the perpendicular axis.

    Examples
    --------
    >>> df = pd.DataFrame({'x': [0, 1, 2, 0, 1, 2], 'y': [0, 0, 0, 10, 10, 10]})
    >>> df_with_lines = assign_line_ids(df, angle_deg=0, spacing=10)
    >>> df_with_lines['line_id'].unique()
    array([0, 1])
    """
    df_out = df.copy()

    x = df[x_col].values
    y = df[y_col].values

    # Compute perpendicular coordinate
    # For a line at angle theta, perpendicular direction is theta + 90
    theta_rad = np.deg2rad(angle_deg)
    perp_rad = theta_rad + np.pi / 2

    # Project onto perpendicular axis
    n_perp = x * np.cos(perp_rad) + y * np.sin(perp_rad)

    # Estimate number of lines
    perp_extent = n_perp.max() - n_perp.min()
    n_clusters = int(np.round(perp_extent / spacing)) + 1
    n_clusters = max(2, n_clusters)

    # 1D KMeans clustering
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=0)
    labels = km.fit_predict(n_perp.reshape(-1, 1))

    # Remap labels to spatial order
    centers = km.cluster_centers_.flatten()
    order = np.argsort(centers)
    remap = {old: new for new, old in enumerate(order)}
    line_id = np.array([remap[label] for label in labels])

    df_out['line_id'] = line_id

    print(f"Assigned {n_clusters} line IDs at angle {angle_deg}° with spacing ~{spacing}")

    return df_out


def leave_lines_out(
    df_xyvtcs: pd.DataFrame,
    holdout_ids: list[int],
    line_angle_deg: float = 0.0,
    line_spacing: float = 50.0,
    x_col: str = 'x',
    y_col: str = 'y',
    value_col: str = 'value',
    set_col: str = 'set',
    df_fl_with_ids: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Hold out flight lines from the MGSIM input for cross-validation.

    Identifies distinct flight lines among observation points, removes the
    specified lines from conditioning (by setting value=NaN and set=0), and
    returns the reduced DataFrame alongside the held-out ground truth.

    Parameters
    ----------
    df_xyvtcs : pd.DataFrame
        Full MGSIM input DataFrame with columns [x, y, value, trend, cluster, set].
        Observation points have set=1 and non-NaN values.
    holdout_ids : list[int]
        Line IDs to hold out (as assigned by assign_line_ids).
    line_angle_deg : float, default=0.0
        Flight line angle in degrees for line identification (0 = E-W).
    line_spacing : float, default=50.0
        Approximate spacing between flight lines (in coordinate units).
    x_col, y_col : str
        Coordinate column names.
    value_col : str
        Value column name.
    set_col : str
        Set indicator column name (1=observation, 0=grid-only).
    df_fl_with_ids : pd.DataFrame, optional
        Pre-computed flight line DataFrame with 'line_id' column
        (e.g., from a previous call to assign_line_ids). If provided,
        line_angle_deg and line_spacing are ignored.

    Returns
    -------
    df_reduced : pd.DataFrame
        Copy of df_xyvtcs with held-out points converted to grid-only
        (value=NaN, set=0). Grid geometry is preserved.
    df_holdout : pd.DataFrame
        The held-out observation points with their original values.
        Includes 'line_id' column.
    df_fl_ids : pd.DataFrame
        All flight line points with assigned line IDs (useful for plotting
        and subsequent calls).

    Notes
    -----
    The held-out points remain in the grid (same row count) but are treated
    as non-observation points by MGSIM. This preserves the grid geometry
    while removing their conditioning influence.

    Examples
    --------
    >>> df_reduced, df_holdout, df_fl_ids = leave_lines_out(
    ...     df_xyvtcs, holdout_ids=[10, 25],
    ...     line_angle_deg=0, line_spacing=50
    ... )
    >>> # Run MGSIM on reduced data
    >>> df_mgsim = mgsim(mg_resols, df_reduced, df_gamma, ...)
    >>> # Compare realizations at held-out locations
    """
    # Identify flight line points
    obs_mask = df_xyvtcs[set_col] == 1

    if df_fl_with_ids is not None:
        df_fl_ids = df_fl_with_ids.copy()
    else:
        # Assign line IDs to observation points only
        df_obs = df_xyvtcs.loc[obs_mask].copy()
        df_fl_ids = assign_line_ids(
            df_obs, angle_deg=line_angle_deg, spacing=line_spacing,
            x_col=x_col, y_col=y_col
        )

    n_lines = df_fl_ids['line_id'].nunique()
    n_obs = len(df_fl_ids)

    # Validate holdout IDs
    available_ids = set(df_fl_ids['line_id'].unique())
    invalid_ids = set(holdout_ids) - available_ids
    if invalid_ids:
        raise ValueError(
            f"Holdout line IDs {invalid_ids} not found. "
            f"Available IDs: 0 to {max(available_ids)}"
        )

    # Extract held-out points
    holdout_mask_fl = df_fl_ids['line_id'].isin(holdout_ids)
    df_holdout = df_fl_ids.loc[holdout_mask_fl].copy()
    holdout_indices = df_holdout.index

    # Build reduced DataFrame
    df_reduced = df_xyvtcs.copy()
    df_reduced.loc[holdout_indices, value_col] = np.nan
    df_reduced.loc[holdout_indices, set_col] = 0

    # Report
    n_holdout = len(df_holdout)
    n_remaining = obs_mask.sum() - n_holdout
    print(f"Leave-lines-out: {n_lines} total lines, holding out {len(holdout_ids)} "
          f"({n_holdout} points)")
    print(f"  Remaining observations: {n_remaining}")
    print(f"  Held-out lines: {holdout_ids}")

    return df_reduced, df_holdout, df_fl_ids


def extract_ensemble_along_line(
    ds: xr.Dataset,
    df_line: pd.DataFrame,
    value_var: str = 'newtrend',
    x_col: str = 'x',
    y_col: str = 'y',
) -> np.ndarray:
    """
    Extract ensemble realization values along a line of points.

    Parameters
    ----------
    ds : xr.Dataset
        Ensemble dataset with dims (realization, point) or similar.
        Must have x, y coordinates and the value variable.
    df_line : pd.DataFrame
        Line points with coordinates to extract.
    value_var : str, default='newtrend'
        Name of the variable in ds containing simulated values.
    x_col, y_col : str
        Coordinate column names in df_line.

    Returns
    -------
    np.ndarray
        Extracted values with shape (n_realizations, n_line_points).

    Notes
    -----
    Uses pandas merge on coordinates for spatial matching.
    Assumes integer or close-to-integer coordinates that match exactly.
    """
    n_realizations = ds.dims.get('realization', len(ds[value_var]))

    # Get coordinates from dataset
    ds_x = ds['x'].values
    ds_y = ds['y'].values

    # Handle different coordinate shapes
    if ds_x.ndim == 2:
        # (realization, point) - take first realization's coords
        ds_x = ds_x[0, :]
        ds_y = ds_y[0, :]
    elif ds_x.ndim == 1 and 'realization' in ds['x'].dims:
        ds_x = ds['x'].isel(realization=0).values
        ds_y = ds['y'].isel(realization=0).values

    line_points = df_line[[x_col, y_col]].values
    n_line_points = len(line_points)

    result = np.zeros((n_realizations, n_line_points))

    for r in range(n_realizations):
        # Get values for this realization
        values = ds[value_var].isel(realization=r).values

        # Create DataFrame for merge
        df_realz = pd.DataFrame({
            x_col: ds_x.ravel(),
            y_col: ds_y.ravel(),
            'val': values.ravel(),
        })

        # Merge with line points
        df_merged = pd.merge(
            df_line[[x_col, y_col]].reset_index(drop=True),
            df_realz,
            on=[x_col, y_col],
            how='left'
        )

        result[r, :] = df_merged['val'].values

    return result


# =============================================================================
# Animation Configuration
# =============================================================================

@dataclass
class ValidationAnimationConfig:
    """Configuration for validation animation."""

    figsize: tuple[float, float] = (16, 14)
    dpi: int = 150
    fps: int = 2
    vmin: float = -1000.0
    vmax: float = 1500.0
    hillshade: bool = True
    realization_alpha: float = 0.15
    realization_color: str = 'steelblue'
    line_highlight_color: str = 'red'
    font_size: int = 12
    title_format: str = "Flight Line {line_id}"
    method_colors: Dict[str, str] = field(default_factory=lambda: {
        'Ground Truth': 'black',
        'Kriging': 'brown',
        'SGSIM Mean': 'green',
        'MGSIM Mean': 'cyan',
    })


# =============================================================================
# Plotting Helpers
# =============================================================================

def plot_rmd_falsification(
    ax: plt.Axes,
    rmd_obs: float,
    rmd_threshold: float,
    rmd_realizations: np.ndarray,
    title: str = 'RMD Falsification',
    quantile: float = 95.0,
) -> None:
    """
    Plot RMD falsification scatter with threshold line.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axes to plot on.
    rmd_obs : float
        RMD value of the observation.
    rmd_threshold : float
        Quantile threshold value.
    rmd_realizations : np.ndarray
        RMD values for each realization.
    title : str
        Plot title.
    quantile : float
        Quantile used for threshold (for label).
    """
    n_realz = len(rmd_realizations)

    # Scatter plot of realization RMDs
    scatter = ax.scatter(
        np.arange(1, n_realz + 1),
        rmd_realizations,
        c=np.abs(rmd_realizations),
        cmap='winter_r',
        s=30,
        vmin=rmd_realizations.min(),
        vmax=rmd_realizations.max(),
        linewidths=0.5,
        edgecolor='k',
        alpha=0.7,
    )

    # Observation marker
    ax.scatter(
        [0],
        [rmd_obs],
        c=[rmd_obs],
        cmap='winter_r',
        marker='D',
        s=100,
        vmin=rmd_realizations.min(),
        vmax=rmd_realizations.max(),
        linewidths=2,
        edgecolor='red',
        zorder=10,
    )

    # Threshold line
    ax.axhline(
        y=rmd_threshold,
        color='red',
        linestyle='--',
        linewidth=2,
        label=f'{quantile}th percentile',
    )

    # Label for observation
    ax.text(
        n_realz * 0.05,
        rmd_obs,
        r'$d_{obs}$',
        color='red',
        fontweight='bold',
        fontsize=10,
        verticalalignment='center',
    )

    ax.set_xlabel('Realization index')
    ax.set_ylabel('Robust Mahalanobis Distance')
    ax.set_xlim(-0.05 * n_realz, 1.05 * n_realz)
    ax.set_title(title)
    ax.legend(loc='upper right', fontsize=9)


def _reshape_to_grid(
    values: np.ndarray,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
) -> tuple[np.ndarray, list[float]]:
    """
    Reshape flat values to 2D grid.

    Returns
    -------
    grid : np.ndarray
        2D array of values.
    extent : list[float]
        [xmin, xmax, ymin, ymax] for imshow.
    """
    x_unique = np.unique(x_coords)
    y_unique = np.unique(y_coords)
    nx, ny = len(x_unique), len(y_unique)

    grid = values.reshape(ny, nx)
    extent = [x_unique.min(), x_unique.max(), y_unique.min(), y_unique.max()]

    return grid, extent


# =============================================================================
# Frame Plotting
# =============================================================================

def plot_validation_frame(
    ground_truth_grid: pd.DataFrame,
    ground_truth_line: pd.DataFrame,
    methods: Dict[str, pd.DataFrame],
    mgsim_ensemble: np.ndarray,
    sgsim_ensemble: np.ndarray,
    flightline_points: pd.DataFrame,
    line_id: int,
    config: ValidationAnimationConfig | None = None,
    fig: plt.Figure | None = None,
) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """
    Plot a single validation frame with 3-row layout.

    Layout:
    - Row 1: Ground truth map with highlighted flight line
    - Row 2: MGSIM realizations + method means | MGSIM RMD
    - Row 3: SGSIM realizations + method means | SGSIM RMD

    Parameters
    ----------
    ground_truth_grid : pd.DataFrame
        Full grid ground truth with columns [x, y, val].
    ground_truth_line : pd.DataFrame
        Ground truth along lines with columns [x, y, val, line_id].
    methods : dict[str, pd.DataFrame]
        Comparison methods: {'Kriging': df, 'SGSIM Mean': df, 'MGSIM Mean': df}.
        Each DataFrame has columns [x, y, val].
    mgsim_ensemble : np.ndarray
        MGSIM realizations for current line, shape (n_realz, n_line_points).
    sgsim_ensemble : np.ndarray
        SGSIM realizations for current line, shape (n_realz, n_line_points).
    flightline_points : pd.DataFrame
        Existing observation locations with columns [x, y].
    line_id : int
        ID of the line to plot.
    config : ValidationAnimationConfig, optional
        Plot configuration.
    fig : plt.Figure, optional
        Existing figure to use.

    Returns
    -------
    fig : plt.Figure
        The figure object.
    axes : dict[str, plt.Axes]
        Dictionary of axes: {'map', 'mgsim_line', 'mgsim_rmd', 'sgsim_line', 'sgsim_rmd'}.
    """
    if config is None:
        config = ValidationAnimationConfig()

    # Get colormap
    cmap = geosoft_cmap_k65()
    norm = Normalize(vmin=config.vmin, vmax=config.vmax)
    ls = LightSource(azdeg=315, altdeg=45)

    # Extract line data
    df_line = ground_truth_line[ground_truth_line['line_id'] == line_id].copy()
    df_line = df_line.sort_values('x').reset_index(drop=True)
    n_line_pts = len(df_line)
    gt_line_vals = df_line['val'].values

    # Create figure if needed
    if fig is None:
        fig = plt.figure(figsize=config.figsize)

    fig.clf()
    gs = GridSpec(3, 2, height_ratios=[1.2, 1, 1], width_ratios=[1.5, 1], figure=fig)

    ax_map = fig.add_subplot(gs[0, :])  # Map spans both columns
    ax_mgsim_line = fig.add_subplot(gs[1, 0])
    ax_mgsim_rmd = fig.add_subplot(gs[1, 1])
    ax_sgsim_line = fig.add_subplot(gs[2, 0])
    ax_sgsim_rmd = fig.add_subplot(gs[2, 1])

    axes = {
        'map': ax_map,
        'mgsim_line': ax_mgsim_line,
        'mgsim_rmd': ax_mgsim_rmd,
        'sgsim_line': ax_sgsim_line,
        'sgsim_rmd': ax_sgsim_rmd,
    }

    # --- Row 1: Map ---
    gt_vals = ground_truth_grid['val'].values
    gt_x = ground_truth_grid['x'].values
    gt_y = ground_truth_grid['y'].values
    gt_grid, extent = _reshape_to_grid(gt_vals, gt_x, gt_y)

    if config.hillshade:
        rgb = ls.shade(gt_grid, cmap=cmap, blend_mode='soft',
                       vmin=config.vmin, vmax=config.vmax)
        ax_map.imshow(rgb, origin='lower', extent=extent, interpolation='nearest')
    else:
        ax_map.imshow(gt_grid, origin='lower', extent=extent, cmap=cmap,
                      norm=norm, interpolation='nearest')

    # Plot existing flight lines (gray)
    ax_map.scatter(
        flightline_points['x'], flightline_points['y'],
        s=3, c='dimgray', alpha=0.5, label='Observations'
    )

    # Highlight current line (red)
    ax_map.scatter(
        df_line['x'], df_line['y'],
        s=30, c=config.line_highlight_color, edgecolors='black',
        linewidths=0.5, zorder=10, label=f'Line {line_id}'
    )

    ax_map.set_title(config.title_format.format(line_id=line_id), fontsize=config.font_size + 2)
    ax_map.set_xlabel('x')
    ax_map.set_ylabel('y')
    ax_map.set_aspect('equal', adjustable='box')
    ax_map.legend(loc='upper right', fontsize=9)

    # Add colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    plt.colorbar(sm, ax=ax_map, shrink=0.6, label='Value')

    # --- Helper function for line plots ---
    def plot_ensemble_comparison(ax, ensemble, ensemble_name, ax_rmd):
        # Plot realizations
        for r in range(len(ensemble)):
            ax.plot(
                ensemble[r],
                color=config.realization_color,
                alpha=config.realization_alpha,
                linewidth=0.8,
            )

        # Plot ground truth
        ax.plot(
            gt_line_vals,
            color=config.method_colors.get('Ground Truth', 'black'),
            linewidth=2.5,
            label='Ground Truth',
        )

        # Plot method means
        for method_name, df_method in methods.items():
            if method_name == 'Ground Truth':
                continue
            # Extract method values along line
            df_method_line = pd.merge(
                df_line[['x', 'y']],
                df_method.rename(columns={'val': 'method_val'}),
                on=['x', 'y'],
                how='left'
            )
            method_vals = df_method_line['method_val'].values
            color = config.method_colors.get(method_name, 'gray')
            ax.plot(method_vals, color=color, linewidth=2, label=method_name)

        ax.set_xlabel('Along-line index')
        ax.set_ylabel('Value')
        ax.set_title(f'{ensemble_name} Realizations vs Methods')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(alpha=0.3, linestyle='--')

        # Compute and plot RMD
        rmd_obs, rmd_thresh, rmd_realz = robust_mahalanobis_distance(
            ensemble, gt_line_vals
        )
        plot_rmd_falsification(
            ax_rmd, rmd_obs, rmd_thresh, rmd_realz,
            title=f'{ensemble_name} RMD (obs={rmd_obs:.2f})'
        )

    # --- Row 2: MGSIM ---
    plot_ensemble_comparison(ax_mgsim_line, mgsim_ensemble, 'MGSIM', ax_mgsim_rmd)

    # --- Row 3: SGSIM ---
    plot_ensemble_comparison(ax_sgsim_line, sgsim_ensemble, 'SGSIM', ax_sgsim_rmd)

    plt.tight_layout()

    return fig, axes


# =============================================================================
# Animation Creation
# =============================================================================

def create_validation_animation(
    ground_truth_grid: pd.DataFrame,
    ground_truth_line: pd.DataFrame,
    methods: Dict[str, pd.DataFrame],
    mgsim_ds: xr.Dataset | np.ndarray,
    sgsim_ds: xr.Dataset | np.ndarray,
    flightline_points: pd.DataFrame,
    output_path: str | Path,
    config: ValidationAnimationConfig | None = None,
    save_frames: bool = True,
    frame_dir: str | Path | None = None,
) -> Path:
    """
    Create validation animation comparing MGSIM and SGSIM.

    Generates an MP4 animation with one frame per flight line, plus optionally
    saves individual PNG frames.

    Parameters
    ----------
    ground_truth_grid : pd.DataFrame
        Full grid ground truth with columns [x, y, val].
    ground_truth_line : pd.DataFrame
        Ground truth along lines with columns [x, y, val, line_id].
    methods : dict[str, pd.DataFrame]
        Comparison methods: {'Kriging': df, 'SGSIM Mean': df, 'MGSIM Mean': df}.
    mgsim_ds : xr.Dataset or np.ndarray
        MGSIM ensemble. If Dataset, extracts values per line.
        If ndarray, shape should be (n_lines, n_realz, n_line_points) or
        pre-extracted per line.
    sgsim_ds : xr.Dataset or np.ndarray
        SGSIM ensemble, same format as mgsim_ds.
    flightline_points : pd.DataFrame
        Existing observation locations with columns [x, y].
    output_path : str or Path
        Output path for MP4 file.
    config : ValidationAnimationConfig, optional
        Animation configuration.
    save_frames : bool, default=True
        Whether to save individual frames as PNG.
    frame_dir : str or Path, optional
        Directory for PNG frames. If None, creates '{output_stem}_frames/'.

    Returns
    -------
    Path
        Path to the created MP4 file.
    """
    if config is None:
        config = ValidationAnimationConfig()

    output_path = Path(output_path)

    # Setup frame directory
    if save_frames:
        if frame_dir is None:
            frame_dir = output_path.parent / f"{output_path.stem}_frames"
        frame_dir = Path(frame_dir)
        frame_dir.mkdir(parents=True, exist_ok=True)
        print(f"Saving frames to: {frame_dir}")

    # Get unique line IDs
    line_ids = sorted(ground_truth_line['line_id'].unique())
    n_lines = len(line_ids)
    print(f"Creating animation for {n_lines} flight lines")

    # Pre-extract ensemble data for each line if xr.Dataset provided
    mgsim_by_line = {}
    sgsim_by_line = {}

    for line_id in line_ids:
        df_line = ground_truth_line[ground_truth_line['line_id'] == line_id]

        if isinstance(mgsim_ds, xr.Dataset):
            mgsim_by_line[line_id] = extract_ensemble_along_line(mgsim_ds, df_line)
        else:
            # Assume already extracted, indexed by line_id or position
            mgsim_by_line[line_id] = mgsim_ds[line_id] if isinstance(mgsim_ds, dict) else mgsim_ds

        if isinstance(sgsim_ds, xr.Dataset):
            sgsim_by_line[line_id] = extract_ensemble_along_line(sgsim_ds, df_line)
        else:
            sgsim_by_line[line_id] = sgsim_ds[line_id] if isinstance(sgsim_ds, dict) else sgsim_ds

    # Create figure
    fig = plt.figure(figsize=config.figsize)

    def update(frame_idx):
        line_id = line_ids[frame_idx]
        print(f"  Frame {frame_idx + 1}/{n_lines}: Line {line_id}")

        plot_validation_frame(
            ground_truth_grid=ground_truth_grid,
            ground_truth_line=ground_truth_line,
            methods=methods,
            mgsim_ensemble=mgsim_by_line[line_id],
            sgsim_ensemble=sgsim_by_line[line_id],
            flightline_points=flightline_points,
            line_id=line_id,
            config=config,
            fig=fig,
        )

        # Save frame
        if save_frames:
            frame_path = frame_dir / f"line_{line_id:03d}.png"
            fig.savefig(frame_path, dpi=config.dpi, bbox_inches='tight')

    # Create animation
    ani = FuncAnimation(
        fig,
        update,
        frames=n_lines,
        interval=1000 // config.fps,
        repeat=True,
    )

    # Save animation
    print(f"Saving animation to: {output_path}")
    writer = FFMpegWriter(fps=config.fps, bitrate=1800)
    ani.save(str(output_path), writer=writer, dpi=config.dpi)

    plt.close(fig)
    print(f"Animation saved: {output_path}")

    return output_path
