# Multigrid Simulation: multigrid-sgsim
This repository houses the multigrid simulation method, a Python package for stochastic interpolation of airborne geophysical data.
This is the software supporting the manuscript:

**Interpolation of large-scale airborne geophysical data with uncertainty quantification**<br>
(Authors: [Joshua H. Rines](https://github.com/jharlanr), [Zhen Yin](https://github.com/sdyinzhen), [Jonas Kloeckner](https://github.com/JonasKloeckner), [Jef Caers](https://profiles.stanford.edu/jef-caers))


published in the journal [Computers \& Geosciences](https://www.sciencedirect.com/science/article/pii/S0098300426000981).

The method presented here combines [automatic segmentation](https://link.springer.com/article/10.1007/s11004-012-9413-6) with multigrid simulation.

## Getting Started
### Installation
Clone the repository and install the package in a clean environment:
```bash
git clone https://github.com/Stanford-Mineral-X/multigrid-sgsim.git
cd multigrid-sgsim
python -m venv venv_mgsim
source venv_mgsim/bin/activate   # (Windows: venv_mgsim\Scripts\activate)
python -m pip install -U pip
pip install -e .
```
Or with conda: `conda env create -f environment.yml && conda activate mgsim`.

### Run demo notebooks
The three notebooks in `demos/` form one seeded, fully reproducible chain — from the raw flight-line data (`demos/data/fl_xyvc.csv`) to the interpolated ensembles:

1. `demo_asm.ipynb` — automatic segmentation (ASM) of the field into stationary subregions; reproduces the shipped cluster assignments exactly.
2. `demo_variogram.ipynb` — trend fitting, normal-score transform, and per-subregion variogram fitting; reproduces the shipped config (`config/config_mgsim_clusteriso_medium.json`) exactly.
3. `demo_interpolation.ipynb` — minimum curvature, kriging, SGSIM, and MGSIM from the same data, ending in an ensemble comparison figure. Set `N_REALIZATIONS` and `SEED` in the top cell; realization *i* is seeded with `SEED + i`, so ensembles are exactly repeatable.

```bash
jupyter notebook demos/demo_interpolation.ipynb   # pip install jupyter if needed
```

## Repo Tree
```
.
├── config
│   └── config_mgsim_clusteriso_medium.json
├── demos
│   ├── data
│   │   ├── config_cluster_medium_data.csv
│   │   ├── fl_xyvc.csv
│   │   └── gt_xyvc.csv
│   ├── demo_asm.ipynb
│   ├── demo_interpolation.ipynb
│   └── demo_variogram.ipynb
├── multigrid_sgsim
│   ├── __init__.py
│   ├── mgsim.py
│   ├── sampling.py
│   ├── segmenting.py
│   ├── trendmaking.py
│   ├── utils.py
│   └── variograms.py
├── environment.yml
├── pyproject.toml
├── README.md
└── tests
    └── test_imports.py
```

## License
This code is released for non-commercial and research purposes. For commercial use, please contact the authors.
