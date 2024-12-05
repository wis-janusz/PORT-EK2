import pathlib
import os
import yaml
import subprocess
import math
import regex
import pandas as pd
from Bio import SeqIO

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
            self.bowtie2_path = config["bowtie2_path"]
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

            self.ref_seq = ".".join(config["ref_seq"].split(".")[:-1])
            if "ref_genes" in config.keys():
                self.ref_genes = config["ref_genes"]
            else:
                self.ref_genes = None

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

        # self.kmer_set = None
        # self.sample_list = None
        # self.sample_group_dict = None
        # self.enriched_groups = None
        # self.avg_cols = [f"{group}_avg" for group in self.sample_groups]
        # self.err_cols = None
        # self.p_cols = None

    def _check_bowtie2_path(self):
        return os.path.exists(f"{self.bowtie2_path}/bowtie2")

    def _check_index_built(self):
        index_files = list(
            pathlib.Path(f"{self.project_dir}/temp/ref_index/").glob(
                f"{self.ref_seq}.*"
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
            f"{self.bowtie2_path}/bowtie2-build",
            "-f",
            f"{self.project_dir}/input/{self.ref_seq}.fasta",
            f"{self.project_dir}/temp/ref_index/{self.ref_seq}",
        ]
        result = subprocess.run(build_cmd, capture_output=True, text=True)
        if verbose == True:
            print(build_cmd)
        if result.returncode != 0:
            raise RuntimeError(result.stderr)
        else:
            if verbose == True:
                print(result.stdout)

    def _bowtie_map(self, verbose: bool = False):
        seed_length = int(math.ceil(self.k / 2))
        map_cmd = [
            f"{self.bowtie2_path}/bowtie2",
            "-a",
            "--norc",
            "--no-hd",
            "-L",
            f"{seed_length}",
            "-x",
            f"{self.project_dir}/temp/ref_index/{self.ref_seq}",
            "-f",
            f"{self.project_dir}/temp/enriched_{self.k}mers.fasta",
            "-S",
            f"{self.project_dir}/temp/enriched_{self.k}mers.sam",
        ]
        if verbose == True:
            print(" ".join(map_cmd))
        result = subprocess.run(map_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(result.stderr)
        else:
            if verbose == True:
                print(result.stdout)

    def run_mapping(self, verbose: bool = False):
        if self._check_bowtie2_path() == False:
            raise FileNotFoundError(
                f"bowtie2 installation not found in {self.bowtie2_path}! Please install bowtie and input its path in portek_config.yaml."
            )
        if self._check_index_built() == False:
            self._bowtie_build_index(verbose)
        self._bowtie_map(verbose)


    def _parse_CIGAR(self, CIGAR_string:str) -> dict:
        n_repeats = regex.findall(r"\d+", CIGAR_string)
        matches = [char for char in CIGAR_string if char.isalpha()]
        CIGAR_list = []
        CIGAR_dict = {m:[] for m in set(matches)}
        for n,m in zip(n_repeats, matches):
            CIGAR_list.extend(int(n)*[m])
        for i, pos in enumerate(CIGAR_list):
            CIGAR_dict[pos].append(i)
        return CIGAR_dict


    def _detect_clip_CIGAR(self, CIGAR_dict:dict) -> bool:
        if "S" in CIGAR_dict.keys() or "H" in CIGAR_dict.keys():
            return True
        else:
            return False
            

    def _read_sam_to_df(self) -> pd.DataFrame:
        in_df = pd.read_csv(f"{self.project_dir}/temp/enriched_{self.k}mers.sam", index_col=0, sep="\t", header = None)
        mappings_df = in_df.loc[:,[3,5,16]]
        mappings_df = mappings_df.rename(columns={3:"ref_pos", 5:"CIGAR",16:"n_mismatch"})
        mappings_df["CIGAR"] = mappings_df["CIGAR"].apply(self._parse_CIGAR)
        mappings_df["n_mismatch"] = mappings_df["n_mismatch"].apply(lambda text: text.split(":")[-1])
        mappings_df.loc[:,["ref_pos", "n_mismatch"]] = mappings_df.loc[:,["ref_pos", "n_mismatch"]].astype(int)
        return mappings_df

    
    def _find_variants(self, mapping:pd. Series):
        if mapping["n_mismatch"] != 0:
            mutations = []
            if ("I" not in mapping["CIGAR"].keys()) and ("D" not in mapping["CIGAR"].keys()):
                q_seq = mapping.name
                t_seq = self.ref_seq[(mapping["ref_pos"]-1):(mapping["ref_pos"]-1+self.k)]
                print(len(q_seq), len(t_seq))
                for i, nuc in enumerate(q_seq):
                    if nuc != t_seq[i]:
                        mutations.append(f"{t_seq[i]}{mapping['ref_pos']+i}{q_seq[i]}")
            return " ".join(mutations)
        else:
            return "WT" 
    

    def analyze_mapping(self):
        mappings_df = self._read_sam_to_df()
        mappings_df = mappings_df.loc[~mappings_df["CIGAR"].apply(self._detect_clip_CIGAR)]
        mappings_df["variants"] = mappings_df.apply(lambda row: self._find_variants(row), axis=1)
        return mappings_df