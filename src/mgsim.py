import numpy as np
import pandas as pd
import gstatsim as gs
from typing import Optional, Union, Dict
from sklearn.preprocessing import QuantileTransformer
from sampling import subsample_dataframe


def mgsim(
    mg_resols,
    df_xyvtcs,
    df_gamma=None,
    vario=None,
    xx: str = 'x',
    yy: str = 'y',
    zz: str = 'residual',
    kk: str = 'cluster',
    num_points: int = 10,
    radius: float = 400,
    # radii: list = [400, 300, 200, 150, 100],
    sgs_or_krige: str = 'sgs',
):

    """
    Perform multigrid sequential Gaussian simulation (MGSGIM) to update the residual field and trend of a dataset.

    Parameters:
    -----------
    mg_resols : list of float
        Sequence of grid spacings for subsampled conditioning data at each multigrid iteration
    df_xyvtcs : pd.DataFrame
        Input dataframe containing the following columns:
        - 'x', 'y' (coordinate locations of each data point)
        - 'value' (measured value at each data point location)
        - 'trend' (initial trend at each point in the grid)
        - 'set' (flag: 1 for observation point location, 0 for non-observation point location)
        - 'cluster' (integer cluster ID or each data point; -1 for locations not to be simulated)
    df_gamma : pd.DataFrame, optional
        Variogram model parameters for cluster-specific simulation (subregions mode).
        One row per cluster. Use this for cluster_sgs. Mutually exclusive with vario.
    vario : list, optional
        Single variogram parameters [azimuth, nugget, major_range, minor_range, sill, vtype]
        for global simulation (no clusters). Use this for okrige_sgs. Mutually exclusive with df_gamma.
    xx, yy, zz, kk : str
        Column names in df_xyvtcs for x, y, residual, and cluster ID respectively
        Defaults: 'x', 'y', 'residual', 'cluster'
    num_points : int
        Number of nearest neighbors to use in SGS step
    radius: float
        Search radius for nearest neighbors in SGS step

    Returns:
    --------
    df_all : pd.DataFrame
        DataFrame with updated 'newtrend' column after multigrid SGSIM; retains all original rows and indices.
        'newtrend' is NaN for rows where 'cluster' < 0 (not simulated).

    Notes:
    ------
    Must provide exactly one of df_gamma or vario:
    - df_gamma: Uses cluster_sgs (cluster-specific variograms)
    - vario: Uses okrige_sgs (single global variogram)
    """
    # Validate inputs: must provide exactly one of df_gamma or vario
    if df_gamma is None and vario is None:
        raise ValueError("Must provide either df_gamma (for subregions) or vario (for global)")
    if df_gamma is not None and vario is not None:
        raise ValueError("Provide only one of df_gamma or vario, not both")

    use_global = vario is not None

    # keep a copy of ALL rows
    df_all = df_xyvtcs.copy()

    # simulate only where cluster >= 0 (preserve original index!)
    sim_mask = df_all[kk] >= 0
    df = df_all.loc[sim_mask].copy()

    # define the locations of simulation (grid and observation points)
    pred_xy_grid = df[['x','y']].values

    # compute initial residual
    df['residual'] = df['value'] - df['trend']

    # initialize a 'newtrend' column
    df['newtrend'] = df['trend'].copy()

    # loop over resolutions
    for (i,mg_resol) in enumerate(mg_resols):
        print(f"MultiGrid iteration {i+1}: Processing resolution {mg_resol}")

        # get sub-dataframe just at observation points ('set' column = 1)
        df_obspts = df[df['set'] == 1].copy()
        
        # mg sample residuals at set resolution
        if i<len(mg_resols)-1:
            df_mgsmpl = subsample_dataframe(df_obspts, column_for_sampling='residual', spacing=mg_resol)
            # print(f" subsmampled to {len(df_mgsmpl)} points, search radius = {1*mg_resol} ")
        else:
            df_mgsmpl = df_obspts.copy()  # last iteration uses all obspts
            # print(f" last iteration, using all {len(df_mgsmpl)} observation points (no subsampling), search radius = {1*mg_resol} ")

        # sequential gaussian simulation of residuals subset
        if sgs_or_krige == 'sgs':
            if use_global:
                # Single variogram for entire field (no clusters) - use okrige_sgs
                mgsgs = gs.Interpolation.okrige_sgs(pred_xy_grid, df_mgsmpl, xx, yy, zz, num_points, vario, radius)
            else:
                # Cluster-specific variograms - use cluster_sgs
                mgsgs = gs.Interpolation.cluster_sgs(pred_xy_grid, df_mgsmpl, xx, yy, zz, kk, num_points, df_gamma, radius)
        elif sgs_or_krige == 'krige':
            # Ordinary kriging (no stochastic component)
            if use_global:
                mgsgs, _ = gs.Interpolation.okrige(pred_xy_grid, df_mgsmpl, xx, yy, zz, num_points, vario, radius)
            else:
                # For kriging with clusters, extract first variogram
                vario_k = df_gamma['Variogram'][0]
                mgsgs, _ = gs.Interpolation.okrige(pred_xy_grid, df_mgsmpl, xx, yy, zz, num_points, vario_k, radius) 

        # update trend on grid and obspts (trend = trend + simulated_residuals)
        df['newtrend'] = df['newtrend'] + mgsgs # NOTE THIS MAYBE WE COULD MAKE MORE ROBUST TO ENSURE THAT WE ARE ADDING THE RIGHT SIMULATED VALUE AT THE RIGHT LOCATION TO THE TREND THERE

        # update residuals on obspts (residuals = data - trend)
        df['residual'] = df['value'] - df['newtrend']

    # # stitch back to full frame (newtrend = NaN where cluster<0)
    # df_all['newtrend'] = np.nan
    # df_all.loc[sim_mask, 'newtrend'] = df['newtrend']

    # stitch back to full frame (newtrend, residual = NaN where cluster<0)
    df_all['newtrend'] = np.nan
    df_all['residual'] = np.nan  # final (post-MG) residual
    df_all.loc[sim_mask, 'newtrend'] = df['newtrend']
    df_all.loc[sim_mask, 'residual'] = df['residual'] 


    return df_all


def mgsim_nst(
    mg_resols,
    df_xyvtcs,
    df_gamma=None,
    vario=None,
    xx: str = 'x',
    yy: str = 'y',
    zz: str = 'Nresidual',
    kk: str = 'cluster',
    num_points: int = 10,
    radius: float = 400,
    sgs_or_krige: str = 'sgs',
    nst_trans: Optional[Union[QuantileTransformer, Dict[int, QuantileTransformer]]] = None,
    clip_nst: bool = True,
    clip_percentile: float = 100.0,
    debug: bool = False,
):
    """
    Perform multigrid sequential Gaussian simulation (MGSGIM) with Normal Score Transform.

    Parameters:
    -----------
    mg_resols : list of float
        Sequence of grid spacings for subsampled conditioning data at each multigrid iteration
    df_xyvtcs : pd.DataFrame
        Input dataframe containing the following columns:
        - 'x', 'y' (coordinate locations of each data point)
        - 'value' (measured value at each data point location)
        - 'trend' (initial trend at each point in the grid)
        - 'set' (flag: 1 for observation point location, 0 for non-observation point location)
        - 'cluster' (integer cluster ID or each data point; -1 for locations not to be simulated)
    df_gamma : pd.DataFrame, optional
        Variogram model parameters for cluster-specific simulation (subregions mode).
        One row per cluster. Use this for cluster_sgs. Mutually exclusive with vario.
    vario : list, optional
        Single variogram parameters [azimuth, nugget, major_range, minor_range, sill, vtype]
        for global simulation (no clusters). Use this for okrige_sgs. Mutually exclusive with df_gamma.
    xx, yy, zz, kk : str
        Column names in df_xyvtcs for x, y, residual, and cluster ID respectively
        Defaults: 'x', 'y', 'Nresidual', 'cluster'
    num_points : int
        Number of nearest neighbors to use in SGS step
    radius : float
        Search radius for nearest neighbors in SGS step
    nst_trans : transformer object
        Fitted sklearn transformer (QuantileTransformer or PowerTransformer)
    clip_nst : bool
        If True (default), clip simulated values in normal space before inverse transforming.
    clip_percentile : float
        Percentile for clipping bounds (default 100.0 = use full range).
        Use e.g. 99.0 to clip to 1st-99th percentile of transformed data.
    debug : bool
        If True, print debugging statistics at each iteration.

    Returns:
    --------
    df_all : pd.DataFrame
        DataFrame with updated 'newtrend' column after multigrid SGSIM; retains all original rows and indices.
        'newtrend' is NaN for rows where 'cluster' < 0 (not simulated).
    """
    # Validate inputs: must provide exactly one of df_gamma or vario
    if df_gamma is None and vario is None:
        raise ValueError("Must provide either df_gamma (for subregions) or vario (for global)")
    if df_gamma is not None and vario is not None:
        raise ValueError("Provide only one of df_gamma or vario, not both")

    use_global = vario is not None

    # keep a copy of ALL rows
    df_all = df_xyvtcs.copy()

    # simulate only where cluster >= 0 (preserve original index!)
    sim_mask = df_all[kk] >= 0
    df = df_all.loc[sim_mask].copy()

    # define the locations of simulation (grid and observation points)
    pred_xy_grid = df[['x','y']].values

    # compute initial residual
    df['residual'] = df['value'] - df['trend']

    # normal score the initial residual and store clip bounds
    data2norm = df['residual'].values.reshape(-1, 1)
    df['Nresidual'] = nst_trans.transform(data2norm)

    # store clipping bounds based on percentile (use observation points only)
    obs_mask = df['set'] == 1
    if clip_nst:
        obs_Nresidual = df.loc[obs_mask, 'Nresidual'].values
        obs_residual = df.loc[obs_mask, 'residual'].values
        if clip_percentile >= 100.0:
            nst_min = obs_Nresidual.min()
            nst_max = obs_Nresidual.max()
            res_min = obs_residual.min()
            res_max = obs_residual.max()
        else:
            lower_pct = (100.0 - clip_percentile) / 2.0
            upper_pct = 100.0 - lower_pct
            nst_min = np.percentile(obs_Nresidual, lower_pct)
            nst_max = np.percentile(obs_Nresidual, upper_pct)
            res_min = np.percentile(obs_residual, lower_pct)
            res_max = np.percentile(obs_residual, upper_pct)

    if debug:
        print("=== INITIAL STATE ===")
        print(f"  Original residuals: min={df['residual'].min():.2f}, max={df['residual'].max():.2f}, "
              f"mean={df['residual'].mean():.2f}, std={df['residual'].std():.2f}")
        print(f"  Nresidual: min={df['Nresidual'].min():.3f}, max={df['Nresidual'].max():.3f}, "
              f"mean={df['Nresidual'].mean():.3f}, std={df['Nresidual'].std():.3f}")
        if clip_nst:
            print(f"  Clipping bounds (normal space): [{nst_min:.3f}, {nst_max:.3f}]")
            print(f"  Clipping bounds (original space): [{res_min:.2f}, {res_max:.2f}]")
            print(f"  (percentile={clip_percentile})")

    # initialize a 'newtrend' column
    df['newtrend'] = df['trend'].copy()

    # loop over resolutions
    for (i, mg_resol) in enumerate(mg_resols):
        print(f"MultiGrid iteration {i+1}: Processing resolution {mg_resol}")

        # get sub-dataframe just at observation points ('set' column = 1)
        df_obspts = df[df['set'] == 1].copy()

        # mg sample residuals at set resolution
        if i < len(mg_resols) - 1:
            df_mgsmpl = subsample_dataframe(df_obspts, column_for_sampling='residual', spacing=mg_resol)
        else:
            df_mgsmpl = df_obspts.copy()  # last iteration uses all obspts

        if debug:
            print(f"  Conditioning data: {len(df_mgsmpl)} points")
            print(f"  Nresidual (conditioning): min={df_mgsmpl['Nresidual'].min():.3f}, "
                  f"max={df_mgsmpl['Nresidual'].max():.3f}, mean={df_mgsmpl['Nresidual'].mean():.3f}")

        # sequential gaussian simulation of residuals subset
        if sgs_or_krige == 'sgs':
            if use_global:
                # Single variogram for entire field (no clusters) - use okrige_sgs
                mgsgs = gs.Interpolation.okrige_sgs(pred_xy_grid, df_mgsmpl, xx, yy, zz, num_points, vario, radius)
            else:
                # Cluster-specific variograms - use cluster_sgs
                mgsgs = gs.Interpolation.cluster_sgs(pred_xy_grid, df_mgsmpl, xx, yy, zz, kk, num_points, df_gamma, radius)
        elif sgs_or_krige == 'krige':
            # Ordinary kriging (no stochastic component)
            if use_global:
                mgsgs, _ = gs.Interpolation.okrige(pred_xy_grid, df_mgsmpl, xx, yy, zz, num_points, vario, radius)
            else:
                vario_k = df_gamma['Variogram'][0]
                mgsgs, _ = gs.Interpolation.okrige(pred_xy_grid, df_mgsmpl, xx, yy, zz, num_points, vario_k, radius)

        # convert to array for processing
        data2inorm = np.asarray(mgsgs).reshape(-1, 1)

        if debug:
            n_below = np.sum(data2inorm < nst_min) if clip_nst else 0
            n_above = np.sum(data2inorm > nst_max) if clip_nst else 0
            print(f"  SGS output (normal space): min={data2inorm.min():.3f}, max={data2inorm.max():.3f}, "
                  f"mean={data2inorm.mean():.3f}, std={data2inorm.std():.3f}")
            if clip_nst:
                print(f"  Values outside clip bounds: {n_below} below, {n_above} above "
                      f"({100*(n_below+n_above)/len(data2inorm):.1f}% total)")

        # clip simulated values to observed range before inverse transform
        if clip_nst:
            data2inorm = np.clip(data2inorm, nst_min, nst_max)

        # inverse transform
        mgs_invtrans = nst_trans.inverse_transform(data2inorm).ravel()

        # also clip in original space to prevent residual drift
        if clip_nst:
            mgs_invtrans = np.clip(mgs_invtrans, res_min, res_max)

        if debug:
            print(f"  After inverse transform (clipped): min={mgs_invtrans.min():.2f}, max={mgs_invtrans.max():.2f}, "
                  f"mean={mgs_invtrans.mean():.2f}, std={mgs_invtrans.std():.2f}")

        # update trend on grid and obspts (trend = trend + simulated_residuals)
        df['newtrend'] = df['newtrend'] + mgs_invtrans

        # update residuals on obspts (residuals = data - trend)
        df['residual'] = df['value'] - df['newtrend']

        # update normalized residuals by normal score transforming the updated residuals
        data2norm = df['residual'].values.reshape(-1, 1)
        df['Nresidual'] = nst_trans.transform(data2norm)

        if debug:
            print(f"  Updated residuals: min={df['residual'].min():.2f}, max={df['residual'].max():.2f}, "
                  f"mean={df['residual'].mean():.2f}, std={df['residual'].std():.2f}")
            print(f"  Updated Nresidual: min={df['Nresidual'].min():.3f}, max={df['Nresidual'].max():.3f}, "
                  f"mean={df['Nresidual'].mean():.3f}, std={df['Nresidual'].std():.3f}")
            print(f"  Newtrend: min={df['newtrend'].min():.2f}, max={df['newtrend'].max():.2f}, "
                  f"mean={df['newtrend'].mean():.2f}")
            print()

    # stitch back to full frame
    df_all['newtrend'] = np.nan
    df_all['residual'] = np.nan
    df_all['Nresidual'] = np.nan
    df_all.loc[sim_mask, 'newtrend'] = df['newtrend'].values
    df_all.loc[sim_mask, 'residual'] = df['residual'].values
    df_all.loc[sim_mask, 'Nresidual'] = df['Nresidual'].values

    return df_all