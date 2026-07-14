# test_imports.py

def test_imports():
    """
    check that required libraries and package modules import successfully
    """

    # core
    import numpy, pandas, xarray, matplotlib, scipy, sklearn, skimage, joblib

    # geostats
    import skgstat, gstatsim

    # local package
    import multigrid_sgsim
    from multigrid_sgsim import mgsim, sampling, trendmaking, utils, variograms
    from multigrid_sgsim.utils import geosoft_cmap_k65
    from multigrid_sgsim.variograms import cluster_variogram
    from multigrid_sgsim.mgsim import mgsim as mgsim_fn, mgsim_nst

    # smoke test
    geosoft_cmap_k65()
    assert callable(cluster_variogram)
    assert callable(mgsim_fn) and callable(mgsim_nst)
