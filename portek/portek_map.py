import pathlib
import os
import shutil
import yaml
import pickle
import subprocess
import math
import regex
import itertools
import operator
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from Bio import SeqIO

import portek


class MappingPipeline:

    def __init__(self, project_dir: str, k):
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(k) != int:
            raise TypeError("k must by an integer!")
        else:
            self.k = k

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

            self.ref_seq_name = ".".join(config["ref_seq"].split(".")[:-1])
            self.ref_seq = str(
                SeqIO.read(
                    f"{project_dir}/input/{config['ref_seq']}", format="fasta"
                ).seq
            )
            if "ref_genes" in config.keys():
                self.ref_genes = config["ref_genes"]
            else:
                self.ref_genes = None

            self.avg_cols = [f"{group}_avg" for group in self.sample_groups]

        except:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )

        self.matrices = {}
        try:
            self.matrices["enriched"] = pd.read_csv(
                f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0
            )

        except:
            raise FileNotFoundError(
                f"No enriched {self.k}-mers table found in {project_dir}output/ ! Please run PORT-EK enriched first!"
            )

        self.mutations = None

    def _check_bowtie2_path(self):
        return shutil.which("bowtie2")

    def _check_index_built(self):
        index_files = list(
            pathlib.Path(f"{self.project_dir}/temp/ref_index/").glob(
                f"{self.ref_seq_name}.*"
            )
        )
        if len(index_files) == 0:
            return False
        else:
            return True

    def _bowtie_build_index(self, verbose: bool = False):
        if os.path.exists(f"{self.project_dir}/temp/ref_index/") == False:
            os.makedirs(f"{self.project_dir}/temp/ref_index")
        build_cmd = [
            f"bowtie2-build",
            "-f",
            f"{self.project_dir}/input/{self.ref_seq_name}.fasta",
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}",
        ]
        result = subprocess.run(build_cmd, capture_output=True, text=True)
        if verbose == True:
            print(build_cmd)
        if result.returncode != 0:
            raise Exception(result.stderr)
        else:
            if verbose == True:
                print(result.stdout)

    def _bowtie_map(self, verbose: bool = False):
        seed_length = int(math.ceil(self.k / 2))
        map_cmd = [
            f"bowtie2",
            "-a",
            "--norc",
            "-L",
            f"{seed_length}",
            "-x",
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}",
            "-f",
            f"{self.project_dir}/temp/enriched_{self.k}mers.fasta",
            "-S",
            f"{self.project_dir}/temp/enriched_{self.k}mers.sam",
        ]
        if verbose == True:
            print(" ".join(map_cmd))
        result = subprocess.run(map_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(result.stderr)
        else:
            if verbose == True:
                print(result.stdout)

    def run_mapping(self, verbose: bool = False):
        if self._check_bowtie2_path() == None:
            raise FileNotFoundError(
                f"bowtie2 not found! Please install bowtie and add it to your PATH!"
            )
        if self._check_index_built() == False:
            self._bowtie_build_index(verbose=verbose)
        self._bowtie_map(verbose=verbose)

    def _parse_CIGAR(self, CIGAR_string: str) -> dict:
        n_repeats = regex.findall(r"\d+", CIGAR_string)
        matches = [char for char in CIGAR_string if char.isalpha()]
        CIGAR_list = []
        for n, m in zip(n_repeats, matches):
            CIGAR_list.extend(int(n) * [m])
        return CIGAR_list

    def _read_sam_to_df(self) -> pd.DataFrame:
        reads = pysam.AlignmentFile(
            f"{self.project_dir}/temp/enriched_{self.k}mers.sam", mode="r"
        )
        read_dict = {
            "kmer": [],
            "flag": [],
            "ref_pos": [],
            "CIGAR": [],
            "n_mismatch": [],
        }
        for read in reads:
            read_dict["kmer"].append(read.query_name)
            read_dict["flag"].append(read.flag)
            read_dict["ref_pos"].append(read.reference_start)
            if read.cigarstring == None:
                read_dict["CIGAR"].append("")
            else:
                read_dict["CIGAR"].append(read.cigarstring)
            if read.has_tag("NM"):
                read_dict["n_mismatch"].append(read.get_tag("NM", False))
            else:
                read_dict["n_mismatch"].append(0)

        mappings_df = pd.DataFrame(read_dict)
        mappings_df["CIGAR"] = mappings_df["CIGAR"].apply(self._parse_CIGAR)
        mappings_df.loc[:, ["flag", "ref_pos", "n_mismatch"]] = mappings_df.loc[
            :, ["flag", "ref_pos", "n_mismatch"]
        ].astype(int)
        mappings_df["group"] = mappings_df["kmer"].apply(
            lambda kmer: self.matrices["enriched"].loc[kmer, "group"]
        )
        return mappings_df

    def _detect_unmapped_CIGAR(self, CIGAR_list: list) -> bool:
        if "S" in CIGAR_list or "H" in CIGAR_list or len(CIGAR_list) == 0:
            return True
        else:
            return False

    def _align_seqs(self, ref_seq, kmer, map_pos, cigar):
        ref_start = map_pos
        aln_len = len(cigar)
        ref_end = ref_start + aln_len
        q_seq = [nuc for nuc in kmer]
        t_seq = [nuc for nuc in ref_seq[ref_start:ref_end]]
        ref_pos = []
        curr_pos = ref_start - 1
        for i, change in enumerate(cigar):
            if change == "D":
                curr_pos += 1
                q_seq.insert(i, "-")
                ref_pos.append(curr_pos)
            elif change == "I":
                t_seq.insert(i, "-")
                ref_pos.append(curr_pos)
            else:
                curr_pos += 1
                ref_pos.append(curr_pos)

        q_seq = q_seq[:aln_len]
        t_seq = t_seq[:aln_len]
        return q_seq, t_seq, ref_pos

    def _join_indels(self, mutations: list):
        subs = []
        ins_pos = []
        del_pos = []

        for mut in mutations:
            if mut[1] == "ins":
                ins_pos.append(mut[0])
            elif mut[1] == "del":
                del_pos.append(mut[0])
            else:
                subs.append(mut)

        ins_pos = set(ins_pos)
        inss = []
        dels = []

        for start_pos in ins_pos:
            muts = [
                mut[2] for mut in mutations if mut[0] == start_pos and mut[1] == "ins"
            ]
            inss.append((start_pos, "ins", "".join(muts)))

        for k, g in itertools.groupby(
            enumerate(sorted(del_pos)), lambda x: x[0] - x[1]
        ):
            del_group_pos = list(map(operator.itemgetter(1), g))
            dels.append((del_group_pos[0], "del", del_group_pos[-1]))

        grouped_muts = subs + inss + dels
        return grouped_muts

    def _find_variants(self, ref_seq, kmer, map_pos, cigar) -> list:
        mutations = []
        q_seq, t_seq, ref_pos = self._align_seqs(ref_seq, kmer, map_pos, cigar)
        if len(q_seq) != len(t_seq) != len(ref_pos):
            raise ValueError(
                f"Improper alingment of k-mer {q_seq} and reference {t_seq}"
            )
        for i in range(len(q_seq)):
            if q_seq[i] != t_seq[i]:
                if q_seq[i] == "-":
                    mutations.append((ref_pos[i] + 1, "del", ref_pos[i] + 1))
                elif t_seq[i] == "-":
                    mutations.append((ref_pos[i] + 1, "ins", q_seq[i]))
                else:
                    mutations.append((ref_pos[i] + 1, t_seq[i], q_seq[i]))

        mutations = self._join_indels(mutations)
        return mutations

    def _mutation_tuple_to_text(self, mutation_as_tuple: tuple) -> str:
        if mutation_as_tuple[1] == "del":
            mutation_as_text = f"{mutation_as_tuple[0]}_{mutation_as_tuple[2]}del"
        elif mutation_as_tuple[1] == "ins":
            mutation_as_text = f"{mutation_as_tuple[0]}_{mutation_as_tuple[0]+1}ins{mutation_as_tuple[2]}"
        else:
            mutation_as_text = (
                f"{mutation_as_tuple[0]}{mutation_as_tuple[1]}>{mutation_as_tuple[2]}"
            )
        return mutation_as_text

    def analyze_mapping(self, verbose: bool = False):
        mappings_df = self._read_sam_to_df()
        mappings_df["mutations"] = "WT"
        mutations_dict = {}
        for row in mappings_df.itertuples():
            if self._detect_unmapped_CIGAR(row.CIGAR) == True:
                mappings_df.loc[row.Index, ["ref_pos", "n_mismatch"]] = 0
                mappings_df.loc[row.Index, "mutations"] = "NA"
                mappings_df.loc[row.Index, "flag"] = 4
            elif row.n_mismatch > 0:
                mutations_as_tuples = self._find_variants(
                    self.ref_seq, row.kmer, row.ref_pos, row.CIGAR
                )
                mappings_df.loc[row.Index, "mutations"] = "; ".join(
                    [self._mutation_tuple_to_text(mut) for mut in mutations_as_tuples]
                )
                for mutation in mutations_as_tuples:
                    mutation_text = self._mutation_tuple_to_text(mutation)
                    if mutation_text not in mutations_dict.keys():
                        mutations_dict[mutation_text] = [row.kmer]
                    else:
                        mutations_dict[mutation_text].append(row.kmer)

        num_kmers = len(self.matrices["enriched"])
        num_primary_mappings = len(mappings_df[mappings_df["flag"] == 0])
        num_secondary_mappings = len(mappings_df[mappings_df["flag"] == 256])
        num_unmapped = len(mappings_df[mappings_df["flag"] == 4])

        if verbose == True:
            print(
                f"\nMapping of {num_kmers} {self.k}-mers resulted in {num_primary_mappings} primary mappings and {num_secondary_mappings} secondary mappings."
            )
            print(f"{num_unmapped} {self.k}-mers couldn't be mapped.")

        self.matrices["mappings"] = mappings_df
        return mutations_dict

    def aggregate_mutations(self, mutations_dict):
        aggregate_dict = {}
        for mutation, kmers in mutations_dict.items():
            aggregate_dict[mutation] = (
                self.matrices["enriched"].loc[kmers, self.avg_cols].max()
            )
        aggregate_df = pd.DataFrame(aggregate_dict).T.fillna(0).sort_index()
        aggregate_df.to_csv(f"{self.project_dir}/temp/mutations.csv")

    def save_mappings_df(self):
        print(f"\nSaving {self.k}-mer mappings.")
        df_to_save = self.matrices["mappings"][
            ["kmer", "ref_pos", "group", "mutations"]
        ].copy()
        df_to_save["ref_pos"] = df_to_save["ref_pos"].apply(
            lambda pos: pos + 1 if pos > 0 else pos
        )
        df_to_save.index.name = "id"
        df_to_save.sort_values("ref_pos", inplace=True)
        df_to_save.to_csv(
            f"{self.project_dir}/output/enriched_{self.k}mers_mappings.csv"
        )


class RefFreePipeline:

    def __init__(self, project_dir: str, k):
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(k) != int:
            raise TypeError("k must by an integer!")
        else:
            self.k = k

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

        except:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )

        self.matrices = {}
        try:
            self.matrices["enriched"] = pd.read_csv(
                f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0
            )

        except:
            raise FileNotFoundError(
                f"No enriched {self.k}-mers table found in {project_dir}output/ ! Please run PORT-EK enriched first!"
            )
        self.group_distros = None
        self.sample_list = None
        self.sample_group_dict = None

    def _get_samples(self, verbose:bool = False):
        sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/indices").glob(
            "*sample_list.pkl"
        )
        sample_list = []
        for filename in sample_list_in_path:
            with open(filename, mode="rb") as in_file:
                partial_list = pickle.load(in_file)
            group = filename.stem.split("_")[0]
            partial_list = [f"{group}_{sample_name}" for sample_name in partial_list]
            sample_list.extend(partial_list)
        sample_group_dict = {
            f"{group}": [
                sample for sample in sample_list if sample.split("_")[0] == f"{group}"
            ]
            for group in self.sample_groups
        }
        self.sample_list = sample_list
        self.sample_group_dict = sample_group_dict
    
    def get_kmer_pos(self, matrix_type: str, verbose: bool = False):
        self._get_samples(verbose)
        print(f"\nGetting {matrix_type} {self.k}-mer position distributions.")
        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{self.k}mer_*_pos_dict.pkl"
        )
        if verbose == True:
            files_counter = 1
            tot_files = len(self.sample_groups)
        group_distros = {}
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                group = filename.stem.split("_")[1]
                temp_dict = pickle.load(in_file)
                temp_dict = {portek.decode_kmer(id, self.k):positions for id, positions in temp_dict.items()}
                group_distros[group] = {}
                for kmer_group in self.matrices["enriched"]["group"].unique():
                    kmers = self.matrices["enriched"][self.matrices["enriched"]["group"] == kmer_group].index
                    filtered_dict = {kmer:Counter(positions) for kmer, positions in temp_dict.items() if kmer in kmers}
                    ratio_dict = {}
                    for kmer, counter in filtered_dict.items():
                        ratio_dict[kmer] = {position:count/len(self.sample_group_dict[group]) for position,count in counter.items()}
                    group_distros[group][kmer_group] = ratio_dict
            if verbose == True:
                print(
                    f"Loaded {self.k}-mer distributions from {files_counter} of {tot_files} groups.",
                    end="\r",
                    flush=True,
                )
                files_counter += 1
        self.group_distros = group_distros

    def calc_distribution_stats(self, verbose:bool=False):
        bout_bout = Counter()
        for fractional_counts in self.group_distros["Bout"]["Bout_enriched"].values():
            bout_bout |= fractional_counts
        bout_bmain = Counter()
        for fractional_counts in self.group_distros["Bout"]["Bmain_enriched"].values():
            bout_bmain |= fractional_counts
        bout_cons = Counter()
        for fractional_counts in self.group_distros["Bout"]["conserved"].values():
            bout_cons |= fractional_counts

        max_x = max([max(bout_bout.keys()),max(bout_bmain.keys()),max(bout_cons.keys())])

        bout_coverage = pd.DataFrame(0,index=range(0,max_x+1), columns=["bout","bmain","conserved"])
        bout_df = pd.DataFrame(bout_bout)
        bout_coverage.update(bout_df)
        print(bout_coverage)