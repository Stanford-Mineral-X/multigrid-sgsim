import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from typing import Iterable


def stratified_sample(
        xyv: np.ndarray,
        strata_labels: np.ndarray,
        choice_per_strata: int,
        *,
        min_var_fraction: float=0.98,
        max_tries: int=50,
        random_state: int | np.random.Generator | None=None,
) -> np.ndarray:
    """
    Draw a stratified sample: uniform count per stratum (capped at available count),
    retrying until sample variance >= min_var_fraction * full variance.

    Parameters:
    -----------
    xyv : (N,3) numpy array
        Columns [x, y, value]
    strata_labels : (N,) numpy array
        Integer strata labels for each point
    choice_per_strata : int
        Number to select in each stratum (capped by stratum size)
    min_var_fraction : float (default=0.98)
        Target fraction of full variance to achieve in subsample
    max_tries : int (default=50)
        Maximum number of attempts to achieve target variance
    random_state : int | Generator | None (default=None)
        Seed or Generator for reproducibility

    Returns:
    --------
    (M,3) numpy array
        Sampled rows [x, y, value]; preserves input dtype
    """
    rng = np.random.default_rng(random_state)

    xyv = np.asarray(xyv)
    labels = np.asarray(strata_labels)
    if xyv.ndim != 2 or xyv.shape[1] < 3:
        raise ValueError("xyv must be a (N,3) array with value in column 3")
    
    data_var = float(np.var(xyv[:,2]))
    if not np.isfinite(data_var):
        raise ValueError("Variance of input values is not finite")
    
    unique = np.unique(labels)
    smpl_var = -np.inf
    best = None

    for _ in range(max_tries):
        picks = []
        for lab in unique:
            mask = (labels == lab)
            n = int(mask.sum())
            if n==0:
                continue
            k = min(choice_per_strata, n)
            idx = rng.choice(np.flatnonzero(mask), size=k, replace=False)
            picks.append(xyv[idx, :3]) # keep columns x,y,v
        if not picks:
            break
        sample = np.vstack(picks)
        smpl_var = float(np.var(sample[:,2])) if sample.shape[0] > 1 else 0.0
        if best is None or smpl_var > float(np.var(best[:,2])):
            best = sample
        if smpl_var >= min_var_fraction * data_var:
            return sample
        
    # fall back to best attempt (don't cast to int, keep original precision)
    return best if best is not None else xyv[:0, :3]



def multigrid_smpl_kdtree_idx(
    indata_xy_val: np.ndarray,
    spacing: float,
    initial_band: float = 0,  # kept for signature parity; unused by design
    *,
    coeff: float = 1.1,
    max_trials: int = 20,
    random_state: int | np.random.Generator | None = None,
) -> np.ndarray:
    """
    Fast multigrid sampler using KD-Tree; returns indices into `indata_xy_val`.
    `spacing` must be in same units as x,y.

    Strategy:
      - Randomly shift a lattice (dx, dy) in [0, spacing]
      - Snap each lattice node to its nearest observed point (KD-Tree)
      - Keep the subset with the largest variance; early-stop if >= coeff * full variance
    """
    rng = np.random.default_rng(random_state)

    xyv_all = np.asarray(indata_xy_val, float)
    mask = np.isfinite(xyv_all).all(axis=1)
    xyv = xyv_all[mask]
    if xyv.size == 0:
        return np.array([], dtype=int)

    X = xyv[:, :2]
    V = xyv[:, 2]
    var_data = float(np.var(V))

    xmin, ymin = X.min(axis=0)
    xmax, ymax = X.max(axis=0)

    tree = cKDTree(X)
    best_idx_local: np.ndarray | None = None
    best_var = -np.inf

    trials = max(1, int(max_trials))
    for _ in range(trials):
        dx = rng.uniform(0, min(spacing, max(1e-9, xmax - xmin)))
        dy = rng.uniform(0, min(spacing, max(1e-9, ymax - ymin)))

        xs = np.arange(xmin + dx, xmax + spacing * 0.51, spacing)
        ys = np.arange(ymin + dy, ymax + spacing * 0.51, spacing)

        if xs.size == 0 or ys.size == 0:
            C = np.array([[0.5 * (xmin + xmax), 0.5 * (ymin + ymax)]], dtype=float)
        else:
            CX, CY = np.meshgrid(xs, ys, indexing="xy")
            C = np.column_stack([CX.ravel(), CY.ravel()])

        _, idx_local = tree.query(C, k=1, workers=-1)
        idx_local = np.unique(idx_local)

        var_rsmpl = float(np.var(V[idx_local])) if idx_local.size > 1 else 0.0
        if var_rsmpl > best_var:
            best_var = var_rsmpl
            best_idx_local = idx_local
        if var_rsmpl >= coeff * var_data:
            break

    if best_idx_local is None:
        return np.flatnonzero(mask)  # conservative fallback

    original_rows = np.flatnonzero(mask)
    return original_rows[best_idx_local]


def subsample_dataframe(
    df: pd.DataFrame,
    *,
    column_for_sampling: str = "value",
    spacing: float = 5,
    initial_band: float = 0,
    coeff: float = 1.1,
    max_trials: int = 20,
    random_state: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """
    Subsample `df` using KD-tree multigrid; returns a copy of the selected rows.

    Expects df to have columns ['x','y', column_for_sampling].
    """
    arr = df[["x", "y", column_for_sampling]].to_numpy()
    sel_idx = multigrid_smpl_kdtree_idx(
        arr,
        spacing=spacing,
        initial_band=initial_band,
        coeff=coeff,
        max_trials=max_trials,
        random_state=random_state,
    )
    return df.iloc[sel_idx].copy().reset_index(drop=True)