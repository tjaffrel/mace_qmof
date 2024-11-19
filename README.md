# QMOFAnalyzer

The `ptable_heatmap_mace.py` python script contains the QMOFAnalyzer which is a Python class designed for analyzing Metal-Organic Framework (MOF) structures and calculating model-specific energy errors. It can compute and visualize the mean absolute errors (MAE) of energies per atom, as well as per element, for various models. 

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

