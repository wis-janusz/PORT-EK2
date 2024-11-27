import os
import pathlib
import yaml
import pickle
import itertools
import multiprocessing
import warnings
import numpy as np
import pandas as pd

# import matplotlib.pyplot as plt
# import seaborn as sns
from collections import defaultdict

# from sklearn.preprocessing import MinMaxScaler
from time import process_time
from Bio import SeqIO

import portek


class KmerFinder:
    """
    KmerFinder:
    """

    def __init__(self, project_dir: str, maxk: int) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(maxk) != int or maxk < 5 or maxk % 2 == 0:
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
                    if len(seq.id.split("|")) > 1:
                        seq.id = seq.id.split("|")[1]
                elif header_format == "ncbi":
                    seq.id = seq.id.split("|")[0][:-1]

                if "/" in seq.id:
                    raise ValueError(
                        "Sequence ids cannot contain '/'. Specify correct header formats in the config file."
                    )
                seq.seq = seq.seq.upper()

            sample_list = [seq.id for seq in seq_list]
            self.seq_lists.append(seq_list)
            samplelist_path = f"{self.project_dir}/input/indices/"
            if os.path.exists(samplelist_path) == False:
                os.makedirs(samplelist_path)
            with open(
                f"{samplelist_path}/{group}_sample_list.pkl", mode="wb"
            ) as out_file:
                pickle.dump(sample_list, out_file, protocol=pickle.HIGHEST_PROTOCOL)

    def _unknown_nuc(self):
        return "X"

    def _find_kmers(self, k: int, group: str, seq_list: int, verbose: bool = False):
        start_time = process_time()
        if verbose == True:
            print(f"Finding all {k}-mers in {group} sequences.", flush=True)
        encoding = defaultdict(self._unknown_nuc)
        encoding["A"] = "00"
        encoding["C"] = "01"
        encoding["G"] = "10"
        encoding["T"] = "11"
        kmer_set = set()

        group_size = len(seq_list)
        avg_dict = {}
        freq_dict = {}

        for seq in seq_list:
            seqid = seq.id
            seq = [encoding[nuc] for nuc in seq.seq]
            kmers_dict = {}
            kmers_pos_dict = {}
            for i in range(0, len(seq) - k + 1):
                kmer = seq[i : i + k]
                if "X" not in kmer:
                    kmer = int("".join(kmer), base=2)
                    kmer_set.add(kmer)
                    if kmer in kmers_dict.keys():
                        kmers_dict[kmer] += 1
                        kmers_pos_dict[kmer].append(i + 1)
                        avg_dict[kmer] += 1 / group_size

                    else:
                        kmers_dict[kmer] = 1
                        kmers_pos_dict[kmer] = [i + 1]
                        if kmer in freq_dict.keys():
                            freq_dict[kmer] += 1 / group_size
                            avg_dict[kmer] += 1 / group_size
                        else:
                            freq_dict[kmer] = 1 / group_size
                            avg_dict[kmer] = 1 / group_size

            with open(
                f"{self.project_dir}/input/indices/{k}mers/{group}_{seqid}_count.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmers_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)
            with open(
                f"{self.project_dir}/input/indices/{k}mers/{group}_{seqid}_pos.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmers_pos_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_set.pkl", mode="wb"
        ) as out_file:
            pickle.dump(kmer_set, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_freq_dict.pkl", mode="wb"
        ) as out_file:
            pickle.dump(freq_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_avg_dict.pkl", mode="wb"
        ) as out_file:
            pickle.dump(avg_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        if verbose == True:
            print(
                f"Done finding all {k}-mers in {group} sequences.",
                flush=True,
            )

        return (k, process_time() - start_time)

    def find_all_kmers(self, n_jobs: int = 4, verbose: bool = False):
        print(
            f"Finding all k-mers of lengths {5} to {self.maxk} in {len(self.sample_groups)} files.",
            flush=True,
        )
        k_to_test = [k for k in range(5, self.maxk + 1, 2)]
        k_to_test.reverse()

        find_kmers_pool_input = []
        for k in k_to_test:
            indices_path = f"{self.project_dir}/input/indices/{k}mers/"
            if os.path.exists(indices_path) == False:
                os.makedirs(indices_path)
            groups = zip(self.sample_groups, self.seq_lists)
            for group, seq_list in groups:
                find_kmers_pool_input.append([k, group, seq_list, verbose])

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            times = pool.starmap(self._find_kmers, find_kmers_pool_input, chunksize=1)

        print("Done finding all k-mers!")

        time_dict = {k: 0 for k in k_to_test}
        for time in times:
            time_dict[time[0]] += time[1]
        return time_dict


class FindOptimalKPipeline:
    """
    FindOptimalKPipeline:
    """

    def __init__(self, project_dir: str, maxk: int, times) -> None:
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(maxk) != int or maxk < 5 or maxk % 2 == 0:
            raise TypeError(
                "Maximum k must by an odd integer not smaller than 5 and not smaller than minimum k!"
            )
        else:
            self.maxk = maxk

        try:
            with open(f"{project_dir}/config.yaml", "r") as config_file:
                config = yaml.safe_load(config_file)
            self.sample_groups = config["sample_groups"]
            self.avg_cols = [f"{group}_avg" for group in self.sample_groups]
            self.freq_cols = [f"{group}_freq" for group in self.sample_groups]
            self.c_cols = [f"{group}_c" for group in self.sample_groups]
            self.f_cols = [f"{group}_f" for group in self.sample_groups]
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
        self.times = times

    def _calc_metrics(self, k: int, verbose: bool = False):

        start_time = process_time()
        if verbose == True:
            print(f"Calculating metrics for {k}-mers.", flush=True)
        kmer_set = set()
        sample_list = []
        kmer_set_in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{k}mer_*_set.pkl"
        )
        sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
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
        kmer_set = list(kmer_set)

        for filename in sample_list_in_path:
            with open(filename, mode="rb") as in_file:
                partial_list = pickle.load(in_file)
            group = filename.stem.split("_")[0]
            partial_list = [f"{group}_{sample_name}" for sample_name in partial_list]
            sample_list.extend(partial_list)

        all_kmer_stat_matrix = pd.DataFrame(
            0.0,
            index=kmer_set,
            columns=[col for col in self.freq_cols] + [col for col in self.avg_cols],
            dtype=np.float64,
        )

        group_sample_dict = {
            f"{group}": [
                sample for sample in sample_list if sample.split("_")[0] == f"{group}"
            ]
            for group in self.sample_groups
        }

        group_len_dict = {
            f"{group}": len(group_sample_dict[group]) for group in self.sample_groups
        }
        tot_samples = len(sample_list)

        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{k}mer_*_avg_dict.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            column_name = "_".join(filename.stem.split("_")[1:-1])
            all_kmer_stat_matrix[column_name] = temp_dict

        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{k}mer_*_freq_dict.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            column_name = "_".join(filename.stem.split("_")[1:-1])
            all_kmer_stat_matrix[column_name] = temp_dict

        all_kmer_stat_matrix = all_kmer_stat_matrix.fillna(0.0)

        for c_col, freq_col, group in zip(
            self.c_cols, self.freq_cols, self.sample_groups
        ):
            all_kmer_stat_matrix[c_col] = round(
                all_kmer_stat_matrix[freq_col] * group_len_dict[group], 0
            )

        all_kmer_stat_matrix["F"] = (
            all_kmer_stat_matrix[self.c_cols].sum(axis=1) / tot_samples
        )

        with np.errstate(divide="ignore"):
            all_kmer_stat_matrix["H"] = np.where(
                all_kmer_stat_matrix["F"] == 1,
                0,
                -(
                    all_kmer_stat_matrix["F"] * np.log2(all_kmer_stat_matrix["F"])
                    + (1 - all_kmer_stat_matrix["F"])
                    * np.log2(1 - all_kmer_stat_matrix["F"])
                ),
            )
        min_F = 2 / tot_samples
        min_H = -(min_F * np.log2(min_F) + (1 - min_F) * np.log2(1 - min_F))
        tot_kmer = len(all_kmer_stat_matrix)
        all_kmer_stat_matrix = all_kmer_stat_matrix.loc[
            all_kmer_stat_matrix["H"] >= min_H
        ]
        sig_kmer = len(all_kmer_stat_matrix)
        if verbose == True:
            print(
                f"{sig_kmer} of {tot_kmer} {k}-mers left after filtering with entropy threshold of {min_H}.",
                flush=True,
            )
        all_kmer_stat_matrix["unique"] = 0
        for kmer in all_kmer_stat_matrix.index:
            if all(
                [
                    all_kmer_stat_matrix.loc[kmer, avg_col]
                    - all_kmer_stat_matrix.loc[kmer, freq_col]
                    == 0
                    for avg_col, freq_col in zip(self.avg_cols, self.freq_cols)
                ]
            ):
                all_kmer_stat_matrix.loc[kmer, "unique"] = 1
        unique_kmer = len(all_kmer_stat_matrix.loc[all_kmer_stat_matrix["unique"] == 1])

        try:
            spec = unique_kmer / sig_kmer
        except(ZeroDivisionError):
            spec = 0
        dt = process_time() - start_time
        mem = float(f"{np.mean(all_kmer_stat_matrix['H']):.2g}")

        if verbose == True:
            print(f"Done calculating metrics for {k}-mers.", flush=True)

        return k, spec, mem, dt, sig_kmer

    def find_optimal_k(self, n_jobs: int = 4, verbose: bool = False):
        print("Finding optimal k.")
        k_to_test = [(k, verbose) for k in range(5, self.maxk + 1, 2)]

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            results = pool.starmap(self._calc_metrics, k_to_test, chunksize=1)

        result_df = pd.DataFrame(
            0,
            index=range(5, self.maxk + 1, 2),
            columns=["spec", "mem", "dt", "spec_rank", "mem_rank", "dt_rank", "score"],
            dtype=float,
        )
        result_df["dt"] = result_df["dt"].astype("object")

        for result in results:
            if result != None:
                if result[4] >= 100:
                    result_df.loc[result[0], "mem"] = result[2]
                    result_df.loc[result[0], "spec"] = round(result[1] * 100)
                    result_df.loc[result[0], "dt"] = result[3]
                else:
                    result_df = result_df.drop(labels=[result[0]])
                    print(
                        f"Fewer than 100 {result[0]}-mers passed the entropy filter, they will be ignored."
                    )

        for k, time in self.times.items():
            if k in result_df.index:
                result_df.loc[k, "dt"] += time
        result_df["dt"] = result_df["dt"].round(3)

        result_df["spec_rank"] = result_df["spec"].rank(method="max")
        result_df["mem_rank"] = result_df["mem"].rank(method="max")
        result_df["dt_rank"] = result_df["dt"].rank(method="max")
        result_df["score"] = (
            result_df["spec_rank"] + result_df["mem_rank"] - result_df["dt_rank"]
        )
        result_df = result_df.sort_values("score", ascending=False)

        with open(
            f"{self.project_dir}/output/k_selection_results.txt", mode="w"
        ) as out_file:
            out_file.write(f"\nHere are the results of optimal k selection:\n")
            for k in result_df.index:
                out_file.write(
                    f"\nk: {k}, Unique kmers: {result_df.loc[k,'spec']}% ({result_df.loc[k,'spec_rank']}), average k-mer entropy: {result_df.loc[k,'mem']} ({result_df.loc[k,'mem_rank']}), CPU time {result_df.loc[k,'dt']} s ({result_df.loc[k,'dt_rank']}), score {result_df.loc[k,'score']}."
                )
            best_k = result_df.idxmax()["score"]
            out_file.write(
                f"\n\nPORT-EK thinks k value of {best_k} is the best for your data!"
                f"\nIt has the {portek.make_ordinal(1+len(result_df)-result_df.loc[best_k, 'spec_rank'])} highest percentage of unique k-mers, {portek.make_ordinal(1+len(result_df)-result_df.loc[best_k, 'mem_rank'])} highest average entorpy and is {portek.make_ordinal(result_df.loc[best_k, 'dt_rank'])} most efficient by CPU time, respectively."
            )
        with open(
            f"{self.project_dir}/output/k_selection_results.txt", mode="r"
        ) as out_file:
            if verbose == True:
                print(out_file.read())
            else:
                lines = out_file.readlines()
                tail = lines[:3] + lines[-2:]
                print(*tail)
