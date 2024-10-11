import os
import pathlib
import yaml
import pickle
import numpy as np
import pandas as pd
from Bio import Align, SeqIO

import portek


class FindOptimalKPipeline:
    """
    FindOptimalKPipeline:
    """

    def __init__(self, project_dir: str, maxk: int) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(maxk) != int:
            raise TypeError("Maximum k must by an integer!")
        else:
            self.maxk = maxk

        try:
            with open(f"{project_dir}/config.yaml", "r") as config_file:
                config = yaml.safe_load(config_file)
            self.sample_groups = config["sample_groups"]
            self.mode = config["mode"]
            if self.mode == "ovr":
                self.goi = config["goi"]
                self.control_groups = self.sample_groups.copy()
                self.control_groups.remove(self.goi)
            else:
                self.goi = None
                self.control_groups = None

        except:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )

    def __repr__(self) -> str:
        pass

    def find_optimal_k(self):
        spec_k = {}
        eff_k = {}
        score_k = {}
        mem_k = {}

        for k in range(5, self.maxk + 2, 2):
            print(f"Calculating metrics for {k}-mers.")
            kmer_set = set()
            sample_list = []
            kmer_set_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(
                f"*{k}mer_set.pkl"
            )
            sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(
                "*sample_list.pkl"
            )

            for filename in kmer_set_in_path:
                with open(filename, mode="rb") as in_file:
                    partial_set = pickle.load(in_file)
                kmer_set.update(partial_set)
            if len(kmer_set) == 0:
                print(f"No {k}-mers found. Skipping.")
                continue

            for filename in sample_list_in_path:
                with open(filename, mode="rb") as in_file:
                    partial_list = pickle.load(in_file)
                group = filename.stem.split("_")[0]
                partial_list = [
                    f"{group}_{sample_name}" for sample_name in partial_list
                ]
                sample_list.extend(partial_list)
            sample_list.sort()

            all_kmer_matrix = pd.DataFrame(
                0, index=list(kmer_set), columns=sample_list, dtype="uint8"
            )
            group_sample_dict = {
                f"{group}": [
                    sample
                    for sample in sample_list
                    if sample.split("_")[0] == f"{group}"
                ]
                for group in self.sample_groups
            }
            in_path = pathlib.Path(f"{self.project_dir}/input/{k}mer_indices").glob(
                "*_count.pkl"
            )

            for filename in in_path:
                with open(filename, mode="rb") as in_file:
                    temp_dict = pickle.load(in_file)
                sample_name = "_".join(filename.stem.split("_")[:-1])
                count_dict = {f"{sample_name}": temp_dict.values()}
                temp_df = pd.DataFrame(
                    count_dict, index=temp_dict.keys(), dtype="uint8"
                )
                all_kmer_matrix.update(temp_df)

            for group in self.sample_groups:
                all_kmer_matrix[f"{group}_avg"] = all_kmer_matrix.loc[
                    :, group_sample_dict[group]
                ].mean(axis=1)

            if 0 in all_kmer_matrix.index:
                all_kmer_matrix = all_kmer_matrix.drop(0)
            mean_count = all_kmer_matrix.loc[:, sample_list].mean(axis=None)
            err_cols = []

            if self.mode == "ovr":
                for j in range(len(self.control_groups)):
                    err_name = f"{self.goi}-{self.control_groups[j]}_err_"
                    err_cols.append(err_name)
                    all_kmer_matrix[err_name] = (
                        all_kmer_matrix[f"{self.goi}_avg"]
                        - all_kmer_matrix[f"{self.control_groups[j]}_avg"]
                    )
            elif self.mode == "ava":
                for j in range(1, len(self.sample_groups)):
                    for i in range(j):
                        err_name = (
                            f"{self.sample_groups[i]}-{self.sample_groups[j]}_err"
                        )
                        err_cols.append(err_name)
                        all_kmer_matrix[err_name] = (
                            all_kmer_matrix[f"{self.sample_groups[i]}_avg"]
                            - all_kmer_matrix[f"{self.sample_groups[j]}_avg"]
                        )

            all_kmer_matrix["RMSE"] = (
                np.sqrt(((all_kmer_matrix[err_cols]) ** 2).mean(axis=1))
                / mean_count
            )

            specificity = all_kmer_matrix["RMSE"].mean()
            efficiency = len(all_kmer_matrix[all_kmer_matrix["RMSE"]>specificity])/len(all_kmer_matrix)
            df_mem = (all_kmer_matrix.memory_usage(index=True, deep=True).sum())/1024/1024/1024
            score = specificity*efficiency
            print(
                f"Done calculating metrics for {k}-mers."
            )
            spec_k[k] = specificity
            eff_k[k] = efficiency
            score_k[k] = score
            mem_k[k] = df_mem

        with open(f"{self.project_dir}/output/k_selection_results", mode="w") as out_file:
            out_file.write("\nHere are the results of optimal k selection:\n")
            for k in spec_k.keys():
                out_file.write(f"\nk: {k}, specificity {round(spec_k[k],4)}, efficiency {round(eff_k[k],4)}, score {round(score_k[k], 4)}, data frame memory size {round(mem_k[k],2)} GB.")

            best_k = max(score_k, key=score_k.get)
            small_score_k = {k:score for k, score in score_k.items() if mem_k[k] < 12}
            if len(small_score_k) > 0:
                best_small_k = max(small_score_k, key=small_score_k.get)
            else:
                best_small_k = None
            if best_k == best_small_k:
                out_file.write(f"\n\nPORT-EK thinks the best k value for your data is {best_k}!")
            else:
                out_file.write(f"\n\nPORT-EK thinks the best k value for your data is {best_k}!")
                out_file.write(f"However, the resulting data frame is larger than 12 GB ({round(mem_k[best_k], 2)} GB)!")
                if best_small_k == None:
                    out_file.write("Unfortunetaly, no k value produces a data frame smaller than 12 GB.")
                else:
                    out_file.write(f"The best k value that produces a data frame smaller than 12 GB is {best_small_k}.")

        with open(f"{self.project_dir}/output/k_selection_results", mode="r") as out_file:
            print(out_file.read())
            

