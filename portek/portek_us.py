import os
import pathlib
import yaml
import pickle
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats, spatial, cluster
from sklearn import preprocessing, decomposition

import portek


class UnsupervisedKmersPipeline:
    """
    UnsupervisedKmersPipeline:
    """

    def __init__(self, project_dir: str, k: int) -> None:
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

            self.ref_rec = config["ref_seq"]
            if "ref_genes" in config.keys():
                self.ref_genes = config["ref_genes"]
            else:
                self.ref_genes = None
            self.freq_cols = [f"{group}_freq" for group in self.sample_groups]
            self.avg_cols = [f"{group}_avg" for group in self.sample_groups]
            self.c_cols = [f"{group}_c" for group in self.sample_groups]
            self.f_cols = [f"{group}_f" for group in self.sample_groups]
        except:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )

        self.kmer_set = None
        self.sample_list = None
        self.sample_group_dict = None
        self.matrices = {}


    def get_kmers(self, verbose: bool = False):

        kmer_set = set()
        sample_list = []
        kmer_set_in_path = f"{self.project_dir}/input/indices/{self.k}mer_set.pkl"
        sample_list_in_path = f"{self.project_dir}/input/indices/sample_list.pkl"

        with open(kmer_set_in_path, mode="rb") as in_file:
            kmer_set = pickle.load(in_file)
        kmer_set = list(kmer_set)

        with open(sample_list_in_path, mode="rb") as in_file:
            sample_list = pickle.load(in_file)

        if len(kmer_set) == 0:
            raise FileNotFoundError(
                f"No {self.k}-mers found in project directory! Make sure you generate them using PORT-EK find_k."
            )
        tot_samples = len(sample_list)
        self.kmer_set = kmer_set
        self.sample_list = sample_list
        print(f"\nImported {len(kmer_set)} kmers and {len(sample_list)} samples.")

        all_kmer_matrix = pd.DataFrame(
            0.0,
            index=kmer_set,
            columns=["F", "avg"],
            dtype=np.float64,
        )

        avg_dict_path = f"{self.project_dir}/input/indices/{self.k}mer_avg_dict.pkl"
        freq_dict_path = f"{self.project_dir}/input/indices/{self.k}mer_freq_dict.pkl"

        with open(avg_dict_path, mode="rb") as in_file:
            avg_dict = pickle.load(in_file)
        all_kmer_matrix["avg"] = avg_dict

        with open(freq_dict_path, mode="rb") as in_file:
            freq_dict = pickle.load(in_file)
        all_kmer_matrix["F"] = freq_dict

        all_kmer_matrix = all_kmer_matrix.fillna(0.0)

        with np.errstate(divide="ignore"):
            all_kmer_matrix["H"] = np.where(
                all_kmer_matrix["F"] == 1,
                0,
                -(
                    all_kmer_matrix["F"] * np.log2(all_kmer_matrix["F"])
                    + (1 - all_kmer_matrix["F"])
                    * np.log2(1 - all_kmer_matrix["F"])
                ),
            )
        min_F = 0.1
        min_H = -(min_F * np.log2(min_F) + (1 - min_F) * np.log2(1 - min_F))
        print(f"{len(all_kmer_matrix[(all_kmer_matrix['H'] < min_H) & (all_kmer_matrix['F'] < 0.9)])} low entropy rare {self.k}-mers rejected.")
        print(f"{len(all_kmer_matrix[all_kmer_matrix['F'] >= 0.9])} conserved {self.k}-mers rejected.")
        common_kmer_matrix = all_kmer_matrix.loc[
            (all_kmer_matrix["H"] >= min_H) & (all_kmer_matrix["F"] < 0.9)
        ]

        print(
            f"{len(common_kmer_matrix)} {self.k}-mers pass entropy filter at H min of {min_H} and are not >90% conserved."
        )

        print(f"\nGetting {self.k}-mer counts.")
        common_count_matrix = pd.DataFrame(
            0,
            index=common_kmer_matrix.index,
            columns=self.sample_list,
            dtype=np.uint8,
        )

        if verbose == True:
            counter = 1
            tot_files = len(self.sample_list)
        in_path = pathlib.Path(f"{self.project_dir}/input/indices/{self.k}mers").glob(
            "*_count.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            sample_name = "_".join(filename.stem.split("_")[:-1])
            count_dict = {f"{sample_name}": temp_dict.values()}
            temp_df = pd.DataFrame(count_dict, index=temp_dict.keys(), dtype="uint8")
            common_count_matrix.update(temp_df)
            if verbose == True:
                print(
                    f"Loaded {self.k}-mers from {counter} of {tot_files} samples.",
                    end="\r",
                    flush=True,
                )
                counter += 1
        common_kmer_matrix = pd.concat([common_count_matrix, common_kmer_matrix], axis=1)
        del common_count_matrix
        print("\nClustering k-mers based on counts and sequences.")
        common_kmer_matrix.index = common_kmer_matrix.index.map(lambda id: portek.decode_kmer(id, self.k))
        common_kmer_matrix["contig"],_ = portek.cluster_kmers(common_kmer_matrix[self.sample_list])
        print(f"Found {len(common_kmer_matrix['contig'].unique())} unique {self.k}-mer contigs.")

        self.matrices["common"] = common_kmer_matrix
        self.matrices["clusters"] = common_kmer_matrix.groupby("contig").first()

    def plot_volcanos(self, matrix_type):
        print(f"\nPlotting and saving volcano plots of enriched {self.k}-mers.")
        for i in range(len(self.err_cols)):
            err = self.err_cols[i]
            group1 = err.split("_")[0].split("-")[0]
            group2 = err.split("_")[0].split("-")[1]
            logp = f"-log10_{self.p_cols[i]}"
            fig, ax = plt.subplots()
            plt.axhline(y=-np.log10(0.01), color="black")
            plt.axvline(x=0.1, color="black", linestyle="--")
            plt.axvline(x=-0.1, color="black", linestyle="--")
            ax.autoscale()
            ax.set_title(f"{group1} vs {group2} volcano plot")
            ax.set_xlabel("Average kmer count change")
            ax.set_ylabel("-log10 of p-value")
            fig.tight_layout()
            sns.scatterplot(
                data=self.matrices[matrix_type],
                x=err,
                y=logp,
                s=10,
                linewidth=0,
                hue="group",
            )
            plt.savefig(
                f"{self.project_dir}/output/{err}_{matrix_type}_{self.k}mers_volcano.svg",
                dpi=600,
                format="svg",
                bbox_inches="tight",
            )

    def get_counts_for_classifier(self, verbose):
        enriched_kmers = (
            self.matrices["enriched"]
            .loc[self.matrices["enriched"]["group"] != "conserved"]
            .index.map(lambda seq: portek.encode_kmer(seq))
        )
        counts_for_classifier = pd.DataFrame(
            0, index=enriched_kmers, columns=self.sample_list, dtype="uint8"
        )
        print(f"\nGetting enriched {self.k}-mer counts.")

        if verbose == True:
            counter = 1
            tot_files = len(self.sample_list)
        in_path = pathlib.Path(f"{self.project_dir}/input/indices/{self.k}mers").glob(
            "*_count.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            sample_name = "_".join(filename.stem.split("_")[:-1])
            count_dict = {f"{sample_name}": temp_dict.values()}
            temp_df = pd.DataFrame(count_dict, index=temp_dict.keys(), dtype="uint8")
            counts_for_classifier.update(temp_df)
            if verbose == True:
                print(
                    f"Loaded {self.k}-mers from {counter} of {tot_files} samples.",
                    end="\r",
                    flush=True,
                )
                counter += 1

        counts_for_classifier.index = counts_for_classifier.index.map(
            lambda id: portek.decode_kmer(id, self.k)
        )
        counts_for_classifier = counts_for_classifier.T
        self.matrices["counts"] = counts_for_classifier
        print(f"\nSaving {self.k}-mer counts for classification.")
        counts_for_classifier["sample_group"] = counts_for_classifier.index.map(
            lambda name: name.split("_")[0]
        )
        counts_for_classifier.to_csv(
            f"{self.project_dir}/output/{self.k}mer_counts_for_classifier.csv",
            index_label="sample_name",
        )


    def plot_PCA(self):
        print(f"\nPlotting and saving PCA plot of enriched {self.k}-mers.")
        pca = decomposition.PCA(2)
        X_PCA = pca.fit_transform(self.matrices["counts"].drop("sample_group", axis=1))
        y_names = self.matrices["counts"]["sample_group"]
        fig, ax = plt.subplots()
        fig.tight_layout()
        ax.set_title("PCA of k-mer counts")
        ax.set_xlabel("Principal component 1")
        ax.set_ylabel("Principal component 2")
        sns.scatterplot(x=X_PCA[:, 0], y=X_PCA[:, 1], hue=y_names, s=20, linewidth=0)
        plt.savefig(
            f"{self.project_dir}/output/{self.k}mer_PCA.svg",
            dpi=600,
            format="svg",
            bbox_inches="tight",
        )


    def save_counts_for_classifier(self):
        print(f"\nSaving {self.k}-mer counts for classification.")
        counts_for_classifier = (
            self.matrices["enriched"]
            .loc[self.matrices["enriched"]["group"] != "conserved", self.sample_list]
            .T
        )
        counts_for_classifier["sample_group"] = counts_for_classifier.index.map(
            lambda name: name.split("_")[0]
        )
        self.matrices["counts"] = counts_for_classifier
        counts_for_classifier.to_csv(
            f"{self.project_dir}/output/{self.k}mer_counts_for_classifier.csv",
            index_label="sample_name",
        )


    def save_matrix(self, matrix_type: str, full: bool = False):
        print(f"\nSaving {matrix_type} {self.k}-mers matrix.")
        if full == True:
            out_filename = f"{self.project_dir}/output/{matrix_type}_{self.k}mers.csv"
            self.matrices[matrix_type].to_csv(out_filename, index_label="kmer")
        else:
            out_filename = (
                f"{self.project_dir}/output/{matrix_type}_{self.k}mers_stats.csv"
            )
            self.matrices[matrix_type].drop(self.sample_list, axis=1).to_csv(
                out_filename, index_label="kmer"
            )
