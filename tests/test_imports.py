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
    from multigrid_sgsim import mgsim, sampling, trendmaking, utils
    from multigrid_sgsim.utils import geosoft_cmap_k65, cluster_variogram

    # smoke test
    geosoft_cmap_k65()
    assert callable(cluster_variogram)
    assert hasattr(mgsim, "__name__")  # trivial sanity