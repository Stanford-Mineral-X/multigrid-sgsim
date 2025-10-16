# Multigrid Simulation
This repository houses the multigrid simulation method, a Python package for stochastic interpolation of airborne geophysical data.

## Getting Started
### Installation
```bash
pip install "git+https://github.com/Stanford-Mineral-X/multigrid-sgsim.git@main#egg=multigrid-sgsim"
```
## Directory Tree
```
├── examples
│   ├── data
│   │   ├── fl_xyvc.csv
│   │   ├── gt_xyvc.csv
│   │   ├── synthetic_flightlines.csv
│   │   └── synthetic_magnetic_groundtruth.csv
│   ├── demo_asm.ipynb
│   └── demo_mgsim.ipynb
├── pyproject.toml
├── README.md
├── src
│   └── multigrid_sgsim
│       ├── __init__.py
│       ├── mgsim.py
│       ├── sampling.py
│       ├── segmenting.py
│       ├── trendmaking.py
│       └── utils.py
└── tests
    └── test_imports.py
```

## License
This code is released for non-commercial and research purposes. For commercial use, please contact the authors.
