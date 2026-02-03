# src/multigrid_sgsim/__init__.py

from .utils import geosoft_cmap_k65
from .mgsim import mgsim, mgsim_nst
from .segmenting import asm_energy, asm_cluster
from .validation import (
    robust_mahalanobis_distance,
    assign_line_ids,
    extract_ensemble_along_line,
    plot_validation_frame,
    create_validation_animation,
    plot_rmd_falsification,
    ValidationAnimationConfig,
)

__all__ = [
    # Utils
    "geosoft_cmap_k65",
    # Core MGSIM
    "mgsim",
    "mgsim_nst",
    # Segmentation
    "asm_energy",
    "asm_cluster",
    # Validation
    "robust_mahalanobis_distance",
    "assign_line_ids",
    "extract_ensemble_along_line",
    "plot_validation_frame",
    "create_validation_animation",
    "plot_rmd_falsification",
    "ValidationAnimationConfig",
]