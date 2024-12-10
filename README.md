# QMOFAnalyzer

The `ptable_heatmap_mace.py` python script contains the QMOFAnalyzer which is a Python class designed for analyzing Metal-Organic Framework (MOF) structures, calculating energy errors and plot various graphics that are included in the MACE-MP-0 paper. It also can extract the outlier structures of the QMOF for which MACE-MP-0

## Features

- **Energy Error Computation**: Calculate per-atom energy differences between various models, such as `pbe_energy` and `macemp0b_large".
- **MAE Analysis**: Compute element-wise MAE values.
- **Data Visualization**: Generate a series of plots, including:
  - Scatter plots of PBE energy vs. model energy.
  - Heatmaps of element-wise energy errors across the periodic table.
  - Histograms of energy errors as a function of different properties like bandgap, magnetic moment, and topology.

## Installation

```bash
pip install pandas tqdm pymatgen matplotlib numpy

