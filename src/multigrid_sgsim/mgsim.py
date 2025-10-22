import numpy as np
import pandas as pd
import gstatsim as gs
from .sampling import subsample_dataframe


def mgsim(
    mg_resols,
    df_xyvtcs,
    df_gamma,
    xx: str = 'x',
    yy: str = 'y',
    zz: str = 'residual',
    kk: str = 'cluster',
    num_points: int = 10,
    radius: float = 40,
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
    df_gamma : pd.DataFrame
        Variogram model parameters for the simulation (region-specific variograms; see gstatsim documentation for more details)
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
    """
    
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
            print(f" subsmampled to {len(df_mgsmpl)} points ")
        else:
            df_mgsmpl = df_obspts.copy()  # last iteration uses all obspts
            print(f" last iteration, using all {len(df_mgsmpl)} observation points (no subsampling) ")

        # sequential gaussian simulation of residuals subset
        mgsgs = gs.Interpolation.cluster_sgs(pred_xy_grid, df_mgsmpl, xx, yy, zz, kk, num_points, df_gamma, radius)

        # update trend on grid and obspts (trend = trend + simulated_residuals)
        df['newtrend'] = df['newtrend'] + mgsgs # NOTE THIS MAYBRE WE COULD MAKE MORE ROBUST TO ENSURE THAT WE ARE ADDING THE RIGHT SIMULATED VALUE AT THE RIGHT LOCATION TO THE TREND THERE

        # update residuals on obspts (residuals = data - trend)
        df['residual'] = df['value'] - df['newtrend']

    # stitch back to full frame (newtrend = NaN where cluster<0)
    df_all['newtrend'] = np.nan
    df_all.loc[sim_mask, 'newtrend'] = df['newtrend']


    return df_all