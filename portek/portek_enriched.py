import os
import pathlib
import yaml
import pickle
import numpy as np
import pandas as pd
from Bio import Align, SeqIO

import portek


class EnrichedKmersPipeline:
    """
    EnrichedKmersPipeline:
    """

    def __init__(
        self, project_dir: str, k: int, c: float, min_rmse: float, rare_m: int = 0
    ) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(k) != int:
            raise TypeError("k must by an integer!")
        else:
            self.k = k

        if type(c) != float and type(c) != int:
            raise TypeError("c must by a number between 0.0 and 1.0!")
        else:
            if c < 0 or c > 1:
                raise ValueError("c must by a number between 0.0 and 1.0!")
            else:
                self.c = c

        if type(min_rmse) == float or type(min_rmse) == int:
            if min_rmse < 0:
                raise ValueError("min_rmse must by a positive number")
            else:
                self.min_rmse = min_rmse
        else:
            raise ValueError("min_rmse must by a positive number")

        try:
            with open(f"{project_dir}/config.yaml", "r") as config_file:
                config = yaml.safe_load(config_file)
            self.sample_groups = config["sample_groups"]
            self.mode = config["mode"]
            if self.mode == "ovr":
                self.goi = config["goi"]
                self.control_groups = self.sample_groups.copy()
                self.control_groups.remove(self.goi)
            elif self.mode == "ava":
                self.goi = None
                self.control_groups = None
            else:
                raise ValueError(
                    "Unrecognized analysis mode, should by ava or ovr. Check your config file!"
                )

            self.ref_rec = config["ref_seq"]
            if "ref_genes" in config.keys():
                self.ref_genes = config["ref_genes"]
            else:
                self.ref_genes = None
            self.freq_cols = [f"{group}_freq" for group in self.sample_groups]
            self.avg_cols = [f"{group}_avg" for group in self.sample_groups]

        except:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )

        self.kmer_set = None
        self.sample_list = None
        self.sample_group_dict = None
        self.enriched_groups = None
        self.matrices = {
            "rare": None,
            "common": None,
            "rare_similar": None,
            "enriched": None
        }

    def __repr__(self) -> str:
        pass

    def get_kmers(self, save_rare:bool=False, verbose:bool=False):
        kmer_set = set()
        sample_list = []
        kmer_set_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(
            f"*_{self.k}mer_set.pkl"
        )
        sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(
            "*sample_list.pkl"
        )

        for filename in kmer_set_in_path:
            with open(filename, mode="rb") as in_file:
                partial_set = pickle.load(in_file)
            kmer_set.update(partial_set)
        if len(kmer_set) == 0:
            raise FileNotFoundError(
                f"No {self.k}-mers found in project directory! Make sure you generate them using generate_kmers.sh."
            )

        for filename in sample_list_in_path:
            with open(filename, mode="rb") as in_file:
                partial_list = pickle.load(in_file)
            group = filename.stem.split("_")[0]
            partial_list = [f"{group}_{sample_name}" for sample_name in partial_list]
            sample_list.extend(partial_list)
        sample_list.sort()

        all_kmer_matrix = pd.DataFrame(
            0, index=list(kmer_set), columns=sample_list, dtype="uint8"
        )
        sample_group_dict = {
            f"{group}": [
                sample for sample in sample_list if sample.split("_")[0] == f"{group}"
            ]
            for group in self.sample_groups
        }

        self.kmer_set = kmer_set
        self.sample_list = sample_list
        self.sample_group_dict = sample_group_dict
        print(f"Imported {len(kmer_set)} kmers and {len(sample_list)} samples.")

        if verbose == True:
            counter = 1
            tot_files = len(sample_list)
        in_path = pathlib.Path(f"{self.project_dir}/input/{self.k}mer_indices").glob(
            "*_count.pkl"
        )

        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            sample_name = "_".join(filename.stem.split("_")[:-1])
            count_dict = {f"{sample_name}": temp_dict.values()}
            temp_df = pd.DataFrame(count_dict, index=temp_dict.keys(), dtype="uint8")
            all_kmer_matrix.update(temp_df)
            if verbose == True:
                print(
                    f"Loaded {self.k}-mers from {counter} of {tot_files} samples",
                    end="\r",
                    flush=True,
                )
                counter += 1

        all_kmer_matrix.index = all_kmer_matrix.index.map(
            lambda id: portek.decode_kmer(id, self.k)
        )

        bin_kmer_matrix = all_kmer_matrix > 0
        for group in self.sample_groups:
            all_kmer_matrix[f"{group}_freq"] = bin_kmer_matrix.loc[
                :, sample_group_dict[group]
            ].mean(axis=1)
            all_kmer_matrix[f"{group}_avg"] = all_kmer_matrix.loc[
                :, sample_group_dict[group]
            ].mean(axis=1)

        if self.k * "A" in all_kmer_matrix.index:
            all_kmer_matrix = all_kmer_matrix.drop(self.k * "A")

        common_kmer_matrix = portek.filter_kmers(
            all_kmer_matrix, freq_cols=self.freq_cols, cons_thr=self.c
        )
        rare_kmer_matrix = all_kmer_matrix.loc[
            all_kmer_matrix.index[~all_kmer_matrix.index.isin(common_kmer_matrix.index)]
        ]

        self.matrices["common"] = common_kmer_matrix

        if save_rare == True:
            self.matrices["rare"] = rare_kmer_matrix

        print(
            f"\n{len(common_kmer_matrix)} common k-mers remaining after filtering at a threshold of {self.c}."
        )

    def calc_kmer_stats(self, matrix_type: str):
        print(f"Identyfying enriched {self.k}-mers.")
        if self.mode == "ava":
            self.matrices[matrix_type]["seq"] = self.matrices[matrix_type].index
            err_cols = []
            p_cols = []
            for j in range(1, len(self.sample_groups)):
                for i in range(j):
                    err_name = f"{self.sample_groups[i]}-{self.sample_groups[j]}_err"
                    p_name = f"{self.sample_groups[i]}-{self.sample_groups[j]}_p-value"
                    err_cols.append(err_name)
                    p_cols.append(p_name)
                    self.matrices[matrix_type][err_name] = (
                        self.matrices[matrix_type][f"{self.sample_groups[i]}_avg"]
                        - self.matrices[matrix_type][f"{self.sample_groups[j]}_avg"]
                    )
                    self.matrices[matrix_type][p_name] = self.matrices[matrix_type][
                        "seq"
                    ].apply(
                        portek.calc_kmer_pvalue,
                        args=(
                            self.sample_group_dict[self.sample_groups[i]],
                            self.sample_group_dict[self.sample_groups[j]],
                            self.matrices[matrix_type],
                        ),
                    )
                    self.matrices[matrix_type][f"-log10_{p_name}"] = -np.log10(
                        self.matrices[matrix_type][p_name]
                    )
            self.matrices[matrix_type]["RMSE"] = np.sqrt(
                ((self.matrices[matrix_type][err_cols]) ** 2).mean(axis=1)
            )
            self.matrices[matrix_type] = self.matrices[matrix_type].sort_values(
                "RMSE", ascending=False
            )
            self.matrices[matrix_type] = self.matrices[matrix_type].drop("seq", axis=1)
            self.matrices[matrix_type]["group"] = self.matrices[matrix_type].apply(
                portek.assign_kmer_group_ava,
                p_cols=p_cols,
                avg_cols=self.avg_cols,
                freq_cols=self.freq_cols,
                axis=1,
            )
            self.matrices[matrix_type]["exclusivity"] = self.matrices[matrix_type].apply(
                portek.check_exclusivity, avg_cols=self.avg_cols, axis=1
            )

        else:
            self.matrices[matrix_type]["seq"] = self.matrices[matrix_type].index
            err_cols = []
            p_cols = []
            for j in range(len(self.control_groups)):
                err_name = f"{self.goi}-{self.control_groups[j]}_err"
                p_name = f"{self.goi}-{self.control_groups[j]}_p-value"
                err_cols.append(err_name)
                p_cols.append(p_name)
                self.matrices[matrix_type][err_name] = (
                    self.matrices[matrix_type][f"{self.goi}_avg"]
                    - self.matrices[matrix_type][f"{self.control_groups[j]}_avg"]
                )
                self.matrices[matrix_type][p_name] = self.matrices[matrix_type]["seq"].apply(
                    portek.calc_kmer_pvalue,
                    args=(
                        self.sample_group_dict[self.goi],
                        self.sample_group_dict[self.control_groups[j]],
                        self.matrices[matrix_type],
                    ),
                )
                self.matrices[matrix_type][f"-log10_{p_name}"] = -np.log10(
                    self.matrices[matrix_type][p_name]
                )
            self.matrices[matrix_type]["RMSE"] = np.sqrt(
                ((self.matrices[matrix_type][err_cols]) ** 2).mean(axis=1)
            )
            self.matrices[matrix_type] = self.matrices[matrix_type].sort_values(
                "RMSE", ascending=False
            )
            self.matrices[matrix_type] = self.matrices[matrix_type].drop("seq", axis=1)
            self.matrices[matrix_type]["group"] = self.matrices[matrix_type].apply(
                portek.assign_kmer_group_ovr,
                goi=self.goi,
                p_cols=p_cols,
                err_cols=err_cols,
                freq_cols=self.freq_cols,
                axis=1,
            )
            self.matrices[matrix_type]["exclusivity"] = self.matrices[matrix_type].apply(
                portek.check_exclusivity, avg_cols=self.avg_cols, axis=1
            )

        self.enriched_groups = [
            name.split("_")[0]
            for name in self.matrices[matrix_type]["group"].value_counts().index
            if "enriched" in name
        ]

    def save_matrix(self, matrix_type: str, full: bool = False):
        if full == True:
            out_filename = f"{self.project_dir}/output/{matrix_type}_{self.k}mers.csv"
            self.matrices[matrix_type].to_csv(out_filename, index_label="kmer")
        else:
            out_filename = f"{self.project_dir}/output/{matrix_type}_{self.k}mers_stats.csv"
            self.matrices[matrix_type].drop(self.sample_list, axis=1).to_csv(out_filename, index_label="kmer")