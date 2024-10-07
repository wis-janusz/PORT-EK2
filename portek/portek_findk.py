import os
import pathlib
import yaml
import pickle
import pandas as pd
from Bio import Align, SeqIO

import portek

class FindOptimalKPipeline: 
    """
    FindOptimalKPipeline:
    """
    def __init__(self, project_dir:str, maxk:int) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")
        
        if type(maxk) != int:
            raise TypeError("Maximum k must by an integer!")
        else:
            self.maxk = maxk

        self.kmer_set = None
        self.sample_list = None
        self.sample_group_dict = None
        self.common_kmer_matrix = None
        self.rare_kmer_matrix = None

        try:
            with open(f"{project_dir}/config.yaml","r") as config_file:
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

            self.avg_cols = [f"{group}_avg" for group in  self.sample_groups]

        except:
            raise FileNotFoundError(f"No config.yaml file found in directory {project_dir}!")

        def __repr__(self) -> str:
            pass

        def find_optimal_k(self) -> dict:
            for k in range(5,maxk+2,2):
                kmer_set = set()
                sample_list = []
                kmer_set_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(f"*{k}mer_set.pkl")
                sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/").glob("*sample_list.pkl")

                for filename in kmer_set_in_path:
                    with open(filename, mode="rb") as in_file:
                        partial_set = pickle.load(in_file)
                    kmer_set.update(partial_set)

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
                group_sample_dict = {f"{group}":[sample for sample in sample_list if sample.split("_")[0] == f"{group}"] for group in self.sample_groups}

                print(f"\nImported {len(kmer_set)} kmers and {len(sample_list)} samples.")

                #Fill matrix with k-mer counts.
                in_path = pathlib.Path(f"{self.project_dir}/input/{k}mer_indices").glob("*_count.pkl")

                for filename in in_path:
                    with open(filename, mode="rb") as in_file:
                        temp_dict = pickle.load(in_file)
                    sample_name = "_".join(filename.stem.split("_")[:-1])
                    count_dict = {f"{sample_name}": temp_dict.values()}
                    temp_df = pd.DataFrame(count_dict, index=temp_dict.keys(), dtype="uint8")
                    all_kmer_matrix.update(temp_df)

                # Construct a temporary binary count matrix, i.e. one that shows if a k-mer is present in sequence, without regards to actual count.
                # Calculate k-mer frequencies and average counts in groups.
                bin_kmer_matrix = all_kmer_matrix > 0
                for group in self.sample_groups:
                    all_kmer_matrix[f"{group}_freq"] = bin_kmer_matrix.loc[:, group_sample_dict[group]].mean(axis=1)
                    all_kmer_matrix[f"{group}_avg"] = all_kmer_matrix.loc[:, group_sample_dict[group]].mean(axis=1)

                # Remove polyA, as its presence and count is mostly dependant on sequencing quality not viral variant.
                if self.k*"A" in all_kmer_matrix.index:
                    all_kmer_matrix = all_kmer_matrix.drop(self.k*"A")