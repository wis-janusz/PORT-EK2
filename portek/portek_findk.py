import os
import pathlib
import yaml
import pickle
import itertools
import multiprocessing
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from sklearn.preprocessing import MinMaxScaler
from time import process_time
from Bio import SeqIO

import portek


class KmerFinder:
    """
    KmerFinder:
    """

    def __init__(self, project_dir: str, mink: int, maxk: int, steps: list) -> None:
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
        if type(steps) != list:
            raise TypeError("--steps must be a list of positive integers!")
        elif any(s < 1 for s in steps):
            raise ValueError("--steps must be a list of positive integers!")
        else:
            self.steps = steps

        if any([step > self.mink for step in self.steps]):
            warnings.warn("One of the step values is larger than minimal k value. K-mers shorter than step value will not give full coverage of sequences!")
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
                    if len(seq.id.split("|"))>1:
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

    def _find_kmers(
        self, k: int, group: str, seq_list: int, s: int = 1, verbose: bool = False
    ):
        start_time = process_time()
        if verbose == True:
            print(
                f"Finding all {k}-mers in {group} sequences with step {s}.", flush=True
            )
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
            for i in range(0, len(seq) - k + 1, s):
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
                f"{self.project_dir}/input/{s}_step/{k}mer_indices/{group}_{seqid}_count.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmers_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)
            # with open(
            #     f"{self.project_dir}/input/{s}_step/{k}mer_indices/{group}_{seqid}_pos.pkl",
            #     mode="wb",
            # ) as out_file:
            #     pickle.dump(kmers_pos_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/{s}_step/{group}_{k}mer_set.pkl", mode="wb"
        ) as out_file:
            pickle.dump(kmer_set, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/{s}_step/{k}mer_{group}_freq_dict.pkl", mode="wb"
        ) as out_file:
            pickle.dump(freq_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/{s}_step/{k}mer_{group}_avg_dict.pkl", mode="wb"
        ) as out_file:
            pickle.dump(avg_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        if verbose == True:
            print(
                f"Done finding all {k}-mers in {group} sequences with step {s}.",
                flush=True,
            )

        return (s, k, process_time()-start_time)

    def find_all_kmers(self, n_jobs: int = 4, verbose: bool = False):
        print(
            f"Finding all k-mers of lengths {self.mink} to {self.maxk} in {len(self.sample_groups)} files.",
            flush=True,
        )
        k_to_test = [k for k in range(self.mink, self.maxk + 1, 2)]

        find_kmers_pool_input = []
        for k in k_to_test:
            for s in [1,5]:
                indices_path = f"{self.project_dir}/input/{s}_step/{k}mer_indices/"
                if os.path.exists(indices_path) == False:
                    os.makedirs(indices_path)
                groups = zip(self.sample_groups, self.seq_lists)
                for group, seq_list in groups:
                    find_kmers_pool_input.append([k, group, seq_list, s, verbose])

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            times = pool.starmap(self._find_kmers, find_kmers_pool_input, chunksize=1)

        print("Done finding all k-mers!")

        time_dict = {sk:0 for sk in itertools.product(*[[1,5], k_to_test])}
        for time in times:
            time_dict[(time[0], time[1])] += time[2]
        return time_dict


class FindOptimalKPipeline:
    """
    FindOptimalKPipeline:
    """

    def __init__(self, project_dir: str, mink: int, maxk: int, steps: list, times) -> None:
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

        if type(steps) != list:
            raise TypeError("--steps must be a list of positive integers!")
        elif any(s < 1 for s in steps):
            raise ValueError("--steps must be a list of positive integers!")
        else:
            self.steps = steps

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


    def _calc_metrics(self, s: int, k: int, verbose: bool = False):
        
        start_time = process_time()
        if verbose == True:
            print(f"Calculating metrics for {k}-mers with step {s}.", flush=True)
        kmer_set = set()
        sample_list = []
        kmer_set_in_path = pathlib.Path(f"{self.project_dir}/input/{s}_step/").glob(
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
                print(f"No {k}-mers with {s} step found. Skipping.", flush=True)
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
            columns=[col for col in self.freq_cols]+[col for col in self.avg_cols],
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

        in_path = pathlib.Path(f"{self.project_dir}/input/{s}_step/").glob(
            f"{k}mer_*_avg_dict.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            column_name = "_".join(filename.stem.split("_")[1:-1])
            all_kmer_stat_matrix[column_name] = temp_dict

        in_path = pathlib.Path(f"{self.project_dir}/input/{s}_step/").glob(
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
        min_F = (min(group_len_dict.values()) / 2) / tot_samples
        min_H = -(min_F * np.log2(min_F) + (1 - min_F) * np.log2(1 - min_F))

        tot_kmer = len(all_kmer_stat_matrix)
        all_kmer_stat_matrix = all_kmer_stat_matrix.loc[
            all_kmer_stat_matrix["H"] >= min_H
        ]
        sig_kmer = len(all_kmer_stat_matrix)

        all_kmer_stat_matrix["unique"] = 0
        for kmer in all_kmer_stat_matrix.index:
            if all([all_kmer_stat_matrix.loc[kmer, avg_col]-all_kmer_stat_matrix.loc[kmer, freq_col] == 0 for avg_col, freq_col in zip(self.avg_cols, self.freq_cols)]):
                all_kmer_stat_matrix.loc[kmer,"unique"] = 1
        cos_distances = []
        err_cols = []
        if self.mode == "ovr":
            for j in range(len(self.control_groups)):
                err_name = f"{self.goi}-{self.control_groups[j]}_err"
                err_cols.append(err_name)
                avg_counts_goi = all_kmer_stat_matrix[
                    f"{self.goi}_avg"
                ]
                avg_counts_j = all_kmer_stat_matrix[
                    f"{self.sample_groups[j]}_avg"
                ]
                all_kmer_stat_matrix[err_name] = avg_counts_goi - avg_counts_j
                cos_distances.append(
                    1
                    - np.dot(avg_counts_goi, avg_counts_j)
                    / (np.linalg.norm(avg_counts_goi) * np.linalg.norm(avg_counts_j))
                )
        elif self.mode == "ava":
            for j in range(1, len(self.sample_groups)):
                for i in range(j):
                    err_name = f"{self.sample_groups[i]}-{self.sample_groups[j]}_err"
                    err_cols.append(err_name)
                    avg_counts_i = all_kmer_stat_matrix[
                        f"{self.sample_groups[i]}_avg"
                    ]
                    avg_counts_j = all_kmer_stat_matrix[
                        f"{self.sample_groups[j]}_avg"
                    ]
                    all_kmer_stat_matrix[err_name] = avg_counts_i - avg_counts_j
                    cos_distances.append(
                        1
                        - np.dot(avg_counts_i, avg_counts_j)
                        / (np.linalg.norm(avg_counts_i) * np.linalg.norm(avg_counts_j))
                    )

        # max_avg_count = all_kmer_stat_matrix[self.avg_cols].max(axis=None)
        # avg_distance = 0
        # for col in err_cols:
        #     err_norm = all_kmer_stat_matrix[col]/max_avg_count
        #     avg_distance += np.sqrt((err_norm**2).sum())

        # avg_distance = avg_distance

        # all_kmer_stat_matrix["RMSE"] = np.sqrt(
        #     ((all_kmer_stat_matrix[err_cols]) ** 2).mean(axis=1)
        # )/max_avg_count
        
        dt = process_time() - start_time
        mem = sig_kmer/tot_kmer
        # mem = round(all_kmer_stat_matrix.memory_usage(index=True, deep=True).sum()/2**20, 3)
        spec = np.mean(cos_distances)

        # fig, ax = plt.subplots()
        # fig.tight_layout()
        # sns.kdeplot(data=all_kmer_stat_matrix,x='H', y="RMSE", fill=True, cmap="mako")
        # ax.axhline(spec, color="orange")
        # ax.axvline(all_kmer_stat_matrix["H"].mean(), color="orange")
        # ax.set_xlabel("H")
        # ax.set_xlim(0, 1)
        # ax.set_ylabel("Normalized RMSE")
        # ax.set_ylim(0, 1)
        # ax.set_title(f"{s} step {k}-mers normalized RMSE vs H distribution")
        # plt.savefig(
        #     f"{self.project_dir}/temp/{s}step_{k}mer_dist.png", format="png", dpi=300
        # )


        # all_kmer_stat_matrix.index = all_kmer_stat_matrix.index.map(
        #     lambda id: portek.decode_kmer(id, k)
        # )
        # all_kmer_stat_matrix.sort_values("RMSE", ascending=False).to_csv(
        #     f"{self.project_dir}/temp/{s}step_{k}mer_matrix.csv"
        # )

        if verbose == True:
            print(f"Done calculating metrics for {k}-mers with step {s}.", flush=True)
        

        return s, k, spec, mem, dt, sig_kmer

    def find_optimal_k(self, n_jobs: int = 4, verbose: bool = False):
        print("Finding optimal s and k.")
        k_to_test = [k for k in range(self.mink, self.maxk + 1, 2)]
        calc_kmers_pool_input = []
        for k in k_to_test:
            for s in [1,5]:
                calc_kmers_pool_input.append([s,k,verbose])

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            results = pool.starmap(self._calc_metrics, calc_kmers_pool_input, chunksize=1)

        result_df = pd.DataFrame(
            0,
            index=pd.MultiIndex.from_tuples([sk[:-1] for sk in calc_kmers_pool_input]),
            columns=["spec", "mem", "dt", "spec_rank", "mem_rank", "dt_rank", "score"],
            dtype=float,
        )
        result_df["dt"] = result_df["dt"].astype("object")

        for result in results:
            if result != None:
                if result[5] >= 100:
                    result_df.loc[(result[0], result[1]), "mem"] = round(
                        result[3], 4
                    )
                    result_df.loc[(result[0], result[1]), "spec"] = round(result[2], 4)
                    result_df.loc[(result[0], result[1]), "dt"] = result[4]
                else:
                    result_df = result_df.drop(labels=[(result[0], result[1])])
                    print(f"Fewer than 100 {result[1]}-mers with step {result[0]} passed the entropy filter, they will be ignored.")

        for sk, time in self.times.items():
            if sk in result_df.index:
                result_df.loc[sk,"dt"] += time
        result_df["dt"] = result_df["dt"].round(3)

        result_df["spec_rank"] = result_df["spec"].rank(ascending=False)
        result_df["mem_rank"] = result_df["mem"].rank(ascending=False)
        result_df["dt_rank"] = result_df["dt"].rank()
        result_df["score"] = len(result_df)+1-result_df["spec_rank"]-result_df["mem_rank"]-result_df["dt_rank"]
        result_df = result_df.sort_values(["spec_rank", "mem_rank", "dt_rank"])

        with open(
            f"{self.project_dir}/output/k_selection_results.txt", mode="w"
        ) as out_file:
            out_file.write(f"\nHere are the results of optimal step and k selection out of {len(result_df)} combinations:\n")
            for s, k in result_df.index:
                out_file.write(
                    f"\nstep: {s}, k: {k}, average cosine distance: {result_df.loc[(s,k),'spec']} ({result_df.loc[(s,k),'spec_rank']}), kmers passing entropy filter: {result_df.loc[(s,k),'mem']*100}% ({result_df.loc[(s,k),'mem_rank']}), CPU time {result_df.loc[(s,k),'dt']} s ({result_df.loc[(s,k),'dt_rank']}), score {result_df.loc[(s,k),'score']}."
                )
            best_sk = result_df.idxmax()["spec"]
            best_compromise_sk = result_df.idxmax()["score"]
            out_file.write(
                f"\n\nStep and k values of {best_sk} give the best separation of groups in your data!"
                f"\nThey are the {portek.make_ordinal(result_df.loc[best_sk, 'dt_rank'])} and {portek.make_ordinal(result_df.loc[best_sk, 'mem_rank'])} most efficient by CPU and memory usage, respectively"
                f"\n\nStep and k values of {best_compromise_sk} are the best compromise of specificity and efficiency!"
                f"\nThey give the {portek.make_ordinal(result_df.loc[best_compromise_sk, 'spec_rank'])} best speration of groups."
                f"\nThey are the {portek.make_ordinal(result_df.loc[best_compromise_sk, 'dt_rank'])} and {portek.make_ordinal(result_df.loc[best_compromise_sk, 'mem_rank'])} most efficient by CPU and memory usage, respectively"
            )
        with open(
            f"{self.project_dir}/output/k_selection_results.txt", mode="r"
        ) as out_file:
            if verbose == True:
                print(out_file.read())
            else:
                lines = out_file.readlines()
                tail = lines[:3]+lines[-6:]
                print(*tail)
