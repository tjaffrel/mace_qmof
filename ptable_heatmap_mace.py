import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
import matplotlib.pyplot as plt
from pymatviz import ptable_heatmap
import numpy as np

class QMOFAnalyzer:
    def __init__(self, qmof_struct_completed, dict_, e0s, models=None):
        self.qmof_struct_completed = qmof_struct_completed
        self.dict_ = dict_
        self.e0s = e0s
        self.models = models
        self.element_mae_data = defaultdict(lambda: defaultdict(list))
        self.element_frequency = defaultdict(int)
        self.element_mae = {}
        self.mae_data = {}
        self.mae_dfs = {}
        self.dict_atomic = defaultdict(dict)

    def calculate_errors(self):
        for qmof_st in tqdm(self.qmof_struct_completed, desc="Processing"):
            qmof_id_st = qmof_st["qmof_id"]
            try:
                extract_st = qmof_st["structure"]
                struct = Structure.from_dict(extract_st)
                atoms = AseAtomsAdaptor.get_atoms(struct)

                element_counts = struct.composition.get_el_amt_dict()
                total_atoms = sum(element_counts.values())

                e0_isolated = 0.0
                for element, count in element_counts.items():
                    e0_isolated += self.e0s[element][0] * count
                    self.element_frequency[element] += count

                energies_per_atom = {}

                self.dict_atomic[qmof_id_st] = {}

                for model in self.models:
                    energy_key = f"{model}_atomic"
                    energy = self.dict_[qmof_id_st].get(model)
                    if energy is not None:
                        if model != "pbe_energy" and model != "macemof":
                            energy += self.dict_[qmof_id_st].get("pbe_energy_vdw", 0.0)
                        energy -= e0_isolated 
                        energy_per_atom = energy / total_atoms
                        self.dict_atomic[qmof_id_st][energy_key]= energy_per_atom
                        energies_per_atom[model] = energy_per_atom

                pbe_energy_per_atom = energies_per_atom.get("pbe_energy")
                if pbe_energy_per_atom is None:
                    continue

                for model, energy_per_atom in energies_per_atom.items():
                    if model != "pbe_energy":
                        total_error = abs(energy_per_atom - pbe_energy_per_atom)
                        for element, count in element_counts.items():
                            fraction = count / total_atoms
                            self.element_mae_data[element][model].append(total_error * fraction)
            except Exception:
                continue

    def compute_mae(self):
        self.element_mae = {
            element: {
                model: (sum(errors[model]) / self.element_frequency[element] * 1000)
                if self.element_frequency[element] else None
                for model in self.models if model != "pbe_energy"
            }
            for element, errors in self.element_mae_data.items()
        }


    def format_dict_atomic(self):
            formatted_data = []

            for qmof_id, values in self.dict_atomic.items():
                pbe_energy_per_atom = self.dict_atomic[qmof_id].get("pbe_energy_atomic")
                for energy_key, energy_value in values.items():
                    if energy_key == "pbe_energy_atomic" or energy_value is None:
                        continue
                    else:
                        formatted_data.append({
                            "qmof_id": qmof_id,
                            "model": energy_key,
                            "mae_energy_per_atom": abs(energy_value - pbe_energy_per_atom),
                            "model_energy_per_atom": energy_value
                        })

            df = pd.DataFrame(formatted_data)
            return df

    def display_dict_atomic(self):
        df = self.format_dict_atomic()
        return df
    
    def filter_outliers(self, qmof_props):
        df = self.display_dict_atomic()
        merged_df = df.set_index('qmof_id').join(qmof_props, on="qmof_id")

        merged_df['pbe_energy_atomic'] = [
            self.dict_atomic[qid]['pbe_energy_atomic'] for qid in merged_df.index
        ]

        outliers_list = []
        inliers_list = []
        dict_std_error = {}
        x = np.array([data.get("pbe_energy_atomic") for data in self.dict_atomic.values()])
        for model in self.models:
            if model != "pbe_energy":
                y = np.array([data.get(f"{model}_atomic") for data in self.dict_atomic.values()])
                valid_indices = [i for i in range(len(x)) if x[i] is not None and y[i] is not None]
                xp = x[valid_indices]
                yp = y[valid_indices]
                dict_std_error[model] = np.std(xp - yp)
                upper_limit_ori = xp + 3 * dict_std_error[model]
                lower_limit_ori = xp - 3 * dict_std_error[model]
                cond = (yp > upper_limit_ori) | (yp < lower_limit_ori)
                model_pd = merged_df[merged_df["model"] == f"{model}_atomic"]
                outliers = model_pd[cond]
                inliers = model_pd[~cond]

                print(f"Outliers for model {model}: {outliers.shape}")
                print(f"Inliers for model {model}: {inliers.shape}")
                
                outliers_list.append(outliers)
                inliers_list.append(inliers)

        merged_df_outliers = pd.concat(outliers_list, axis=0)
        merged_df_inliers = pd.concat(inliers_list, axis=0)
        
        return merged_df_outliers, merged_df_inliers


    def scatter_plot(self, map_title_name=None):
        from scipy.stats import gaussian_kde
        x = np.array([data.get("pbe_energy_atomic", None) for data in self.dict_atomic.values()])

        for model in self.models:
            if model != "pbe_energy":
                y = np.array([data.get(f"{model}_atomic") for data in self.dict_atomic.values()])
                valid_indices = [i for i in range(len(x)) if x[i] is not None and y[i] is not None]

                xp = np.array([x[i] for i in valid_indices])
                y = np.array([y[i] for i in valid_indices])
                if len(xp) != len(y):
                    continue  

                error = xp - y
                mean_error = np.mean(error)
                std_error = np.std(error)
                print(mean_error, std_error)
                upper_limit_ori = xp + 3 * std_error
                lower_limit_ori = xp - 3 * std_error

                x_min, x_max = np.min(xp), np.max(xp)
                x_range = np.linspace(x_min - 0.3*abs(x_max - x_min), 
                                    x_max + 0.3*abs(x_max - x_min), 100)
                upper_limit = x_range  + 3*std_error
                lower_limit = x_range  - 3*std_error
                outside_count = np.sum((y > upper_limit_ori) | (y < lower_limit_ori))
                total_count = len(y)
                print(f"Model: {model}")
                print(f"Total Points: {total_count}, Points Outside ±3 Std Dev: {outside_count} ({outside_count / total_count:.2%})")

                mae = np.mean(np.abs(xp - y))
                rmse = np.sqrt(np.mean((xp - y) ** 2))
                r2_score = np.corrcoef(xp, y)[0, 1] ** 2

                xy = np.vstack([xp, y])
                kde = gaussian_kde(xy)(xy)
                scatter = plt.scatter(xp, y, c=kde, cmap='viridis', s=10)
                plt.colorbar(scatter, label='# MOFs')
                plt.plot(x_range, x_range, 'r--')
                plt.plot(x_range, upper_limit, color='cyan', linestyle='--')
                plt.plot(x_range, lower_limit, color='cyan', linestyle='--')
                plt.fill_between(
                    x_range, 
                    lower_limit, 
                    upper_limit, 
                    color='lightblue', 
                    alpha=0.3
                )
                plt.xlim(-7, -3)
                plt.ylim(-7, -3)
                plt.xlabel('PBE Energy [eV/atom]')
                if map_title_name is not None:
                    plt.ylabel(f'{map_title_name[model]} Energy [eV/atom]')
                else:
                    plt.ylabel(f'{model} Energy [eV/atom]')
                plt.title(f'QMOF ({int(len(y))} MOF structures)')
                label_text = f'MAE = {mae:.3f} eV\nRMSE = {rmse:.3f} eV\n$r^2$ = {r2_score:.3f}\nOut. ±3 Std = {outside_count} ({outside_count / total_count:.2%})'
                plt.annotate(label_text, xy=(0.05, 0.95), xycoords='axes fraction',
                             ha='left', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.6))
                plt.savefig(f"{model}_scatter_test.png")
                plt.close()

    def generate_heatmap_data(self, dict_diff_models=None):
    
        self.mae_data = {}

        for model in self.models:
            if model != "pbe_energy":
                self.mae_data[model] = {
                    element: round(maes[model], 2)
                    for element, maes in self.element_mae.items() if maes[model] is not None
                }

        if dict_diff_models:
            for diff_model, models_pair in dict_diff_models.items():
                model1, model2 = models_pair
                self.mae_data[diff_model] = {
                    element: round(
                        self.element_mae[element][model1] - self.element_mae[element][model2], 2
                    )
                    for element in self.element_mae
                    if self.element_mae[element].get(model1) is not None and self.element_mae[element].get(model2) is not None
                }

        self.mae_dfs = {
            model: pd.DataFrame(data.items(), columns=["element", model]).set_index("element")
            for model, data in self.mae_data.items()
        }
    

    def plot_heatmap(self, dict_diff_models=None, map_title_name=None):
        from matplotlib.colors import Normalize
        for model in self.mae_dfs:
            if dict_diff_models and model in dict_diff_models:
                cbar_title_ = r"$\Delta E_{MAE}$ $[meV/atom]$"
                max_ = round(self.mae_dfs[model].values.max(), 2)
                min_ = round(self.mae_dfs[model].values.min(), 2)
                cbar_range_ = (-2, 2)
                colormap = 'coolwarm'
            else:
                cbar_title_ = "Energy MAE [meV/atom]"
                max_ = round(self.mae_dfs[model].values.max(), 2)
                cbar_range_ = (0, 4.5)
                colormap = 'viridis'

            norm = Normalize(vmin=cbar_range_[0], vmax=cbar_range_[1])

            fig = ptable_heatmap(
                self.mae_dfs[model],
                log=False,
                anno_kwargs={"fontsize": 8},
                cbar_title=cbar_title_,
                # cbar_kwargs=cbar_kwargs_,
                cbar_range=cbar_range_,
                colormap=colormap,
                return_type="figure",
                cbar_kwargs={"norm": norm}
            )

            title = map_title_name.get(model, model)
            fig.suptitle(
                f"Element-wise MAE {title} for QMOF Structures",
                fontsize=16,
                fontweight="bold"
            )
            fig.savefig(f"{model}_mae_test.png")
            plt.show()



