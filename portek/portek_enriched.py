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
    def __init__(self, project_dir:str, k:int, c:float, min_rmse:float) -> None:
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
            self.ref_rec = config["ref_seq"]
            if "ref_genes" in config.keys():
                self.ref_genes = config["ref_genes"]
            else:
                self.ref_genes = None
            self.freq_cols = [f"{group}_freq" for group in self.sample_groups]
            self.avg_cols = [f"{group}_avg" for group in  self.sample_groups]
            self.aligner = Align.PairwiseAligner(
                scoring="megablast",
                mode="local"
            )
        except:
            raise FileNotFoundError(f"No config.yaml file found in directory {project_dir}!")


        def __repr__(self) -> str:
            pass


        def get_kmers(self, save_rare = False):
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

            self.kmer_set = kmer_set
            self.sample_list = sample_list
            self.group_sample_dict = group_sample_dict
            print(f"\nImported {len(kmer_set)} kmers and {len(sample_list)} samples.")

            #Fill matrix with k-mer counts.
            counter = 1
            tot_files = len(sample_list)
            in_path = pathlib.Path(f"{self.project_dir}/input/{k}mer_indices").glob("*_count.pkl")

            for filename in in_path:
                with open(filename, mode="rb") as in_file:
                    temp_dict = pickle.load(in_file)
                sample_name = "_".join(filename.stem.split("_")[:-1])
                count_dict = {f"{sample_name}": temp_dict.values()}
                temp_df = pd.DataFrame(count_dict, index=temp_dict.keys(), dtype="uint8")
                all_kmer_matrix.update(temp_df)
                print(
                    f"{counter} of {tot_files} indices done.",
                    end="\r",
                    flush=True,
                )
                counter += 1

            # Decode k-mer sequences
            all_kmer_matrix.index = all_kmer_matrix.index.map(lambda id: portek.decode_kmer(id,k))

            # Construct a temporary binary count matrix, i.e. one that shows if a k-mer is present in sequence, without regards to actual count.
            # Calculate k-mer frequencies and average counts in groups.
            bin_kmer_matrix = all_kmer_matrix > 0
            for group in self.sample_groups:
                all_kmer_matrix[f"{group}_freq"] = bin_kmer_matrix.loc[:, group_sample_dict[group]].mean(axis=1)
                all_kmer_matrix[f"{group}_avg"] = all_kmer_matrix.loc[:, group_sample_dict[group]].mean(axis=1)

            # Remove polyA, as its presence and count is mostly dependant on sequencing quality not viral variant.
            if self.k*"A" in all_kmer_matrix.index:
                all_kmer_matrix = all_kmer_matrix.drop(self.k*"A")

            # Apply rarity filter.
            common_kmer_matrix = portek.filter_kmers(
                all_kmer_matrix, freq_cols=self.freq_cols, cons_thr=self.c
            )
            rare_kmer_matrix = all_kmer_matrix.loc[all_kmer_matrix.index[~all_kmer_matrix.index.isin(common_kmer_matrix.index)]]

            # Save matrices to pipeline.
            self.common_kmer_matrix = common_kmer_matrix
            if save_rare == True:
                self.rare_kmer_matrix = rare_kmer_matrix

            print(
                f"\n{len(common_kmer_matrix)} common k-mers remaining after filtering at a threshold of {self.c}."
            )