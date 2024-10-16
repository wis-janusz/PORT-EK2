import os
import pathlib
import yaml
import pickle
import multiprocessing
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio import SeqIO


class KmerFinder:
    """
    KmerFinder:
    """

    def __init__(self, project_dir: str, mink: int, maxk: int) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(mink) != int or mink < 5 or mink % 2 == 0 or mink > maxk:
            raise TypeError(
                "Minimum k must by an odd integer not smaller than 5 and not bigger than maximum k!"
            )
        else:
            self.mink = mink

        if type(maxk) != int or maxk < 5 or maxk % 2 == 0 or maxk < mink:
            raise TypeError(
                "Maximum k must by an odd integer not smaller than 5 and not smaller than minimum k!"
            )
        else:
            self.maxk = maxk

        try:
            with open(f"{project_dir}/config.yaml", "r") as config_file:
                config = yaml.safe_load(config_file)
            self.sample_groups = config["sample_groups"]
            if len(config["header_format"]) != 0:
                self.input_fastas = zip(
                    config["input_files"],
                    config["header_format"],
                    config["sample_groups"],
                )
            else:
                self.input_fastas = zip(
                    config["input_files"],
                    ["plain"] * len(config["input_files"]),
                    config["sample_groups"],
                )

        except:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )

        self.seq_lists = []

        for in_file, header_format, group in self.input_fastas:
            in_path = f"{self.project_dir}/input/{in_file}"
            seq_list = list(SeqIO.parse(in_path, format="fasta"))
            for seq in seq_list:

                if header_format == "gisaid":
                    seq.id = seq.id.split("|")[1]
                elif header_format == "ncbi":
                    seq.id = seq.id.split("|")[0][:-1]

                if "/" in seq.id:
                    raise ValueError(
                        "Sequence ids cannot contain '/'. If using data from GISAID please use '--header_format gisaid' option."
                    )
                seq.seq = seq.seq.upper()

            sample_list = [seq.id for seq in seq_list]
            self.seq_lists.append(seq_list)
            with open(
                f"{self.project_dir}/input/{group}_sample_list.pkl", mode="wb"
            ) as out_file:
                pickle.dump(sample_list, out_file, protocol=pickle.HIGHEST_PROTOCOL)

    def _unknown_nuc(self):
        return "X"

    def _find_kmers(self, k: int, group: str, seq_list: int, verbose: bool = False):
        if verbose == True:
            print(f"Finding all {k}-mers in {group} sequences.", flush=True)
        encoding = defaultdict(self._unknown_nuc)
        encoding["A"] = "00"
        encoding["C"] = "01"
        encoding["G"] = "10"
        encoding["T"] = "11"
        kmer_set = set()

        for seq in seq_list:
            seqid = seq.id
            seq = [encoding[nuc] for nuc in seq.seq]
            kmers_dict = {}
            kmers_pos_dict = {}
            for i in range(len(seq) - k):
                kmer = seq[i : i + k]
                if "X" not in kmer:
                    kmer = int("".join(kmer), base=2)
                    kmer_set.add(kmer)
                    if kmer in kmers_dict.keys():
                        kmers_dict[kmer] += 1
                        kmers_pos_dict[kmer].append(i + 1)
                    else:
                        kmers_dict[kmer] = 1
                        kmers_pos_dict[kmer] = [i + 1]

            with open(
                f"{self.project_dir}/input/{k}mer_indices/{group}_{seqid}_count.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmers_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)
            with open(
                f"{self.project_dir}/input/{k}mer_indices/{group}_{seqid}_pos.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmers_pos_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input//{group}_{k}mer_set.pkl", mode="wb"
        ) as out_file:
            pickle.dump(kmer_set, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        if verbose == True:
            print(f"Done finding all {k}-mers in {group} sequences.")

    def find_all_kmers(self, n_jobs: int = 4, verbose: bool = False):
        print(
            f"Finding all k-mers of lengths {self.mink} to {self.maxk} in {len(self.sample_groups)} files.",
            flush=True,
        )
        k_to_test = [k for k in range(self.mink, self.maxk + 1, 2)]
        find_kmers_pool_input = []
        for k in k_to_test:
            indices_path = f"{self.project_dir}/input/{k}mer_indices/"
            if os.path.exists(indices_path) == False:
                os.mkdir(indices_path)
            groups = zip(self.sample_groups, self.seq_lists)
            for group, seq_list in groups:
                find_kmers_pool_input.append([k, group, seq_list, verbose])

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            pool.starmap(self._find_kmers, find_kmers_pool_input, chunksize=1)
        print("Done finding all k-mers!")


class FindOptimalKPipeline:
    """
    FindOptimalKPipeline:
    """

    def __init__(self, project_dir: str, mink: int, maxk: int) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(mink) != int or mink < 5 or mink % 2 == 0 or mink > maxk:
            raise TypeError(
                "Minimum k must by an odd integer not smaller than 5 and not bigger than maximum k!"
            )
        else:
            self.mink = mink

        if type(maxk) != int or maxk < 5 or maxk % 2 == 0 or maxk < mink:
            raise TypeError(
                "Maximum k must by an odd integer not smaller than 5 and not smaller than minimum k!"
            )
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

    def _calc_metrics(self, k: int, verbose: bool = False):
        if verbose == True:
            print(f"Calculating metrics for {k}-mers.", flush=True)
        kmer_set = set()
        sample_list = []
        kmer_set_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(
            f"*_{k}mer_set.pkl"
        )
        sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/").glob(
            "*sample_list.pkl"
        )

        for filename in kmer_set_in_path:
            with open(filename, mode="rb") as in_file:
                partial_set = pickle.load(in_file)
            kmer_set.update(partial_set)
        if len(kmer_set) == 0:
            if verbose == True:
                print(f"No {k}-mers found. Skipping.", flush=True)
            return None

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
        group_sample_dict = {
            f"{group}": [
                sample for sample in sample_list if sample.split("_")[0] == f"{group}"
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
            temp_df = pd.DataFrame(count_dict, index=temp_dict.keys(), dtype="uint8")
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
                    err_name = f"{self.sample_groups[i]}-{self.sample_groups[j]}_err"
                    err_cols.append(err_name)
                    all_kmer_matrix[err_name] = (
                        all_kmer_matrix[f"{self.sample_groups[i]}_avg"]
                        - all_kmer_matrix[f"{self.sample_groups[j]}_avg"]
                    )

        all_kmer_matrix["RMSE"] = (
            np.sqrt(((all_kmer_matrix[err_cols]) ** 2).mean(axis=1)) / mean_count
        )

        specificity = all_kmer_matrix["RMSE"].mean()
        efficiency = len(all_kmer_matrix[all_kmer_matrix["RMSE"] > specificity]) / len(
            all_kmer_matrix
        )
        df_mem = (
            (all_kmer_matrix.memory_usage(index=True, deep=True).sum())
            / (2**30)
        )
        score = specificity * efficiency
        if verbose == True:
            print(f"Done calculating metrics for {k}-mers.", flush=True)

        return k, specificity, efficiency, score, df_mem

    def find_optimal_k(self, n_jobs: int = 4, verbose:bool = False):
        print("Finding optimal k.")
        spec_k = {}
        eff_k = {}
        score_k = {}
        mem_k = {}
        k_to_test = [(k,verbose) for k in range(self.mink, self.maxk + 1, 2)]
        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            results = pool.starmap(self._calc_metrics, k_to_test, chunksize=1)

        for result in results:
            if result != None:
                spec_k[result[0]] = result[1]
                eff_k[result[0]] = result[2]
                score_k[result[0]] = result[3]
                mem_k[result[0]] = result[4]

        with open(
            f"{self.project_dir}/output/k_selection_results", mode="w"
        ) as out_file:
            out_file.write("\nHere are the results of optimal k selection:\n")
            for k in spec_k.keys():
                out_file.write(
                    f"\nk: {k}, specificity {round(spec_k[k],4)}, efficiency {round(eff_k[k],4)}, score {round(score_k[k], 4)}, data frame memory size {round(mem_k[k],2)} GB."
                )

            best_k = max(score_k, key=score_k.get)
            small_score_k = {k: score for k, score in score_k.items() if mem_k[k] < 12}
            if len(small_score_k) > 0:
                best_small_k = max(small_score_k, key=small_score_k.get)
            else:
                best_small_k = None
            if best_k == best_small_k:
                out_file.write(
                    f"\n\nPORT-EK thinks the best k value for your data is {best_k}!"
                )
            else:
                out_file.write(
                    f"\n\nPORT-EK thinks the best k value for your data is {best_k}!"
                )
                out_file.write(
                    f"\nHowever, the resulting data frame is larger than 12 GB ({round(mem_k[best_k], 2)} GB)!"
                )
                if best_small_k == None:
                    out_file.write(
                        "\nUnfortunetaly, no k value produces a data frame smaller than 12 GB."
                    )
                else:
                    out_file.write(
                        f"\nThe best k value that produces a data frame smaller than 12 GB is {best_small_k}."
                    )

        with open(
            f"{self.project_dir}/output/k_selection_results", mode="r"
        ) as out_file:
            print(out_file.read())
