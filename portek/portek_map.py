import pathlib
import os
import yaml
import subprocess
import pandas as pd
from Bio import SeqIO


class MappingPipeline():

    def __init__(self, project_dir:str, k):
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(k) != int:
            raise TypeError("k must by an integer!")
        else:
            self.k = k
        
        with open("portek_config.yaml", "r") as global_config_file:
            global_config = yaml.safe_load(global_config_file)
        self.bowtie2_path = global_config["bowtie2_path"]

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

            self.ref_seq = config["ref_seq"].split(".")[-1]
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
            self.matrices["enriched"] = pd.read_csv(f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0)

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
        

    def run_mapping(self):
        if pathlib.Path(f"{self.project_dir}/temp/*").match(f"{self.ref_seq}") == False:
            print("no index")
            build_cmd = [f"{self.bowtie2_path}bowtie"]
        cmd = [f"{self.bowtie2_path}bowtie2"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.returncode)
        