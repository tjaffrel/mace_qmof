# QMOFAnalyzer

The `ptable_heatmap_mace.py` python script of the QMOFAnalyzer, this Python class is designed for analyzing MOF structures, calculating energy errors and plot various graphics that are included in the MACE-MP-0 paper. It also can extract the outlier and inlier structures for which MACE-MP-0 on QMOF. 

## Figures
- **Figures 1**: Figure 1a-c; Comparison between the PBE energy per atom with MACE-MP-0b predictions for the different model sizes; (a) large, (b) medium, and (c) small; on the entire QMOF database, 20,375 structures. The blue zone represents the range within 3 standard deviations of the mean absolute error distribution. Structures outside this zone, is only about 2% of the dataset, indicated in the legend. <br />
The mean absolute error (MAE) and root mean square error (RMSE) difference between the large and small model are 10 meV/atom. <br />
The smaller model could thus be used for high-throughput computation for large systems, the ratio size of the model and performance being in favor of the small model. 
- **Figure 2**: Figure 2a-c; Element-wise MAE for the predictions of the MACE-MP-0b models; (a) large, (b) medium, and (c) small; compared to PBE energies on the QMOF database. This Figure is similar to Figure 4b of the article with an improved color scale (difficult to analyze in the original article) in addition of putting the given values in the bottom for a better comparison. I also added a script for better reproducibility. Overall the error results between the three models are very similar. <br />
Figure 2d shows the difference in element-wise MAE between the small and large models. Blue indicate elements better predicted by the small model, respectively red for the large model. The large model exhibit better, or almost similar, performance for almost all elements, except for Technetium, which is rarely used in MOF chemistry.
- **Figure 3**: Figure 3a-c;MAE per atom as a function of atomic density for the MACE-MP-0b models; (a) large, (b) medium, and (c) small; compared to PBE energies. Figure also similar to Figure 43b of your paper but with improved visualization compared to Figure 43b in and enabling comparison among models. <br />
Figure 3d-f; focus on the inlier structures within the blue zone from Figure 1a-c. <br />
Figure 3g-i; focus on the outlier structures outside the blue zone. 
- **Figure 4**: Figure 4a-c; Histograms of the oxidation states of metal secondary building units (SBUs) for 50 MOF in the outliers of Figure 1a-c, for the (a) large, (b) medium, and (c) small models, calculated using the element-based parameters of O'Keeffe et al. <br />
The models exhibit similar distributions of metal SBUs and oxidation states, with a high density for Gd(3+) and rare metals such as V, Mo, and Co, consistent with Figure 2 in the article. The errors for Cu(+), Cu(2+), and Zn(2+) are mainly due to the high number of MOF with such metal oxidation states in MOF chemistry.
- **Figure 5**: Figure 5a-b; 
 Histograms of the averaged MAE per atom across crystal systems for the MOF outliers found in Figure 1a-c, categorized by model size: (a) large, (b) medium, and (c) small.<br />
These figures show that the highest outliers have predominantly orthorhombic, triclinic, and tetragonal crystal cells. The error on the large model appears more uniform across crystal systems, compared to the small and medium ones, suggesting that the errors may mainly only come from the chemistry and the PBE data. 

## Functions

- **Energy Error Computation**: Calculate per-atom energy differences between various models, such as `pbe_energy` and `macemp0b_large".
- **MAE Analysis**: Compute element-wise MAE values.
- **Data Visualization**: Generate a series of plots, including:
  - Scatter plots of PBE energy vs. model energy.
  - Heatmaps of element-wise energy errors across the periodic table.
  - Histograms of energy errors as a function of different properties like crystral systems and band gap.

## Installation

```bash
pip install pandas tqdm pymatgen matplotlib numpy

