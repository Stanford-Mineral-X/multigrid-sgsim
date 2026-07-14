import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skgstat as skg
from typing import List, Optional


def cluster_variogram(df_fl, value_or_residual, idx_cluster, maxlag, n_lags, model,
                      azimuth: Optional[float] = None, bandwidth: Optional[float] = None):
    """
    Compute empirical variogram for a specific cluster, optionally directional.

    Parameters
    ----------
    df_fl : pd.DataFrame
        Input dataframe with coordinates and values
    value_or_residual : str
        Column name for values to compute variogram on
    idx_cluster : int
        Cluster index to filter data
    maxlag : float
        Maximum lag distance
    n_lags : int
        Number of lag bins
    model : str
        Variogram model type. Must be one of 'spherical', 'exponential', or 'gaussian'.
        Note: gstatsim only supports these three model types for SGS.
    azimuth : float, optional
        Direction for directional variogram in degrees (0 = East, 90 = North).
        If None, computes omnidirectional variogram.
    bandwidth : float, optional
        Angular tolerance bandwidth in degrees for directional variogram.
        Only used if azimuth is specified. Default is 15 degrees.

    Returns
    -------
    skg.Variogram
        Fitted variogram object
    """
    df = df_fl[df_fl['cluster'] == idx_cluster]
    coords = df[['x', 'y']].values
    values_for_variogram = df[value_or_residual]

    if azimuth is not None:
        # Directional variogram
        bw = bandwidth if bandwidth is not None else 15.0
        Vgram = skg.DirectionalVariogram(
            coordinates=coords,
            values=values_for_variogram,
            azimuth=azimuth,
            bandwidth=bw,
            model=model,
            bin_func="even",
            n_lags=n_lags,
            maxlag=maxlag,
            normalize=False
        )
    else:
        # Omnidirectional variogram
        Vgram = skg.Variogram(
            coordinates=coords,
            values=values_for_variogram,
            model=model,
            bin_func="even",
            n_lags=n_lags,
            maxlag=maxlag,
            normalize=False
        )

    return Vgram


def build_variogram_dataframe(
    variograms: List[skg.Variogram],
    azimuths: Optional[List[float]] = None,
    minor_ranges: Optional[List[float]] = None,
    vtype: Optional[str] = None,
) -> pd.DataFrame:
    """
    Build variogram parameter dataframe for gstatsim from fitted variogram objects.

    Parameters
    ----------
    variograms : list of skg.Variogram
        List of fitted variogram objects (one per cluster)
    azimuths : float or list of float, optional
        Azimuth angle(s) in degrees for anisotropy direction (0 = East).
        Can be a single value (same for all clusters) or a list (one per cluster).
        If None, defaults to 0 for all clusters.
    minor_ranges : list of float, optional
        Minor ranges for anisotropic variograms (one per cluster).
        If None, uses isotropic (minor_range = major_range).
    vtype : str, optional
        Variogram type string for gstatsim. Must be one of 'Spherical', 'Exponential',
        or 'Gaussian' (case-sensitive). If None, infers from first variogram's model.

    Returns
    -------
    pd.DataFrame
        DataFrame with 'Variogram' column containing [azimuth, nugget, major_range, minor_range, sill, vtype]
    """
    # Infer variogram type from model if not specified
    if vtype is None:
        model_name = variograms[0].model.__name__
        # Capitalize first letter for gstatsim
        vtype = model_name.capitalize()

    # Handle azimuths - can be single value or list
    if azimuths is None:
        azimuths = [0.0] * len(variograms)
    elif isinstance(azimuths, (int, float)):
        azimuths = [float(azimuths)] * len(variograms)

    gams = []
    for i, V in enumerate(variograms):
        # V.parameters = [range, sill, nugget] for spherical/exponential/gaussian
        major_range, sill, nugget = V.parameters

        # Use specified minor range or default to isotropic
        if minor_ranges is not None:
            minor_range = minor_ranges[i]
        else:
            minor_range = major_range

        # Get azimuth for this cluster
        az = azimuths[i]

        # gstatsim format: [azimuth, nugget, major_range, minor_range, sill, variogram_type]
        gam = [az, nugget, major_range, minor_range, sill, vtype]
        gams.append(gam)

    return pd.DataFrame({'Variogram': gams})


def plot_directional_variograms(
    df_fl: pd.DataFrame,
    value_or_residual: str,
    idx_cluster: int,
    maxlag: float,
    n_lags: int,
    model: str,
    azimuths: List[float] = [0, 45, 90, 135],
    bandwidth: float = 22.5,
    figsize: tuple = (12, 8),
) -> List[skg.DirectionalVariogram]:
    """
    Plot directional variograms at multiple azimuths for anisotropy analysis.

    Parameters
    ----------
    df_fl : pd.DataFrame
        Input dataframe with coordinates and values
    value_or_residual : str
        Column name for values to compute variogram on
    idx_cluster : int
        Cluster index to filter data
    maxlag : float
        Maximum lag distance
    n_lags : int
        Number of lag bins
    model : str
        Variogram model type. Must be one of 'spherical', 'exponential', or 'gaussian'.
        Note: gstatsim only supports these three model types for SGS.
    azimuths : list of float
        List of azimuth angles in degrees to compute variograms for
    bandwidth : float
        Angular tolerance bandwidth in degrees
    figsize : tuple
        Figure size for plot

    Returns
    -------
    list of skg.DirectionalVariogram
        Fitted directional variogram objects
    """
    n_dirs = len(azimuths)
    ncols = min(2, n_dirs)
    nrows = (n_dirs + ncols - 1) // ncols

    fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
    if n_dirs == 1:
        axs = [axs]
    else:
        axs = axs.flatten()

    variograms = []
    for i, az in enumerate(azimuths):
        V = cluster_variogram(
            df_fl, value_or_residual, idx_cluster,
            maxlag, n_lags, model,
            azimuth=az, bandwidth=bandwidth
        )
        variograms.append(V)

        # Plot empirical
        axs[i].plot(V.bins, V.experimental, 'o', color='blue', label='Empirical')

        # Plot fitted model
        x_model = np.linspace(0, V.bins[-1], 100)
        y_model = V.transform(x_model)
        axs[i].plot(x_model, y_model, '-', color='red', label=f'Fitted ({model})')

        axs[i].set_xlabel('Lag Distance (h)')
        axs[i].set_ylabel('γ(h)')
        axs[i].set_title(f'Cluster {idx_cluster} - Azimuth {az}°')
        axs[i].legend()
        axs[i].grid(True)

    # Hide unused subplots
    for j in range(i + 1, len(axs)):
        axs[j].set_visible(False)

    plt.tight_layout()
    plt.show()

    return variograms
