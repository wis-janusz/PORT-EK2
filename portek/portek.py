import os
import pathlib
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import regex
import yaml
from matplotlib import pyplot, collections, colormaps, patches, colors
from datetime import datetime
from scipy import stats
from Bio import Align


def encode_kmer(kmer_seq: str) -> int:
    encoding_dict = {"A": "00", "C": "01", "G": "10", "T": "11"}
    kmer_bin_string = [encoding_dict[nuc] for nuc in kmer_seq]
    id = int("".join(kmer_bin_string), base=2)
    return id


def decode_kmer(id: int, k) -> str:
    decoding = {"00": "A", "01": "C", "10": "G", "11": "T"}
    kmer_bin_string = bin(id)[2:].rjust(2 * k, "0")
    kmer_bin_string = [
        kmer_bin_string[i : i + 2] for i in range(0, len(kmer_bin_string), 2)
    ]
    kmer_seq = "".join([decoding[bits] for bits in kmer_bin_string])
    return kmer_seq


def filter_kmers(kmer_df: pd.DataFrame, freq_cols: list, cons_thr=0.01) -> pd.DataFrame:
    out_kmer_df = kmer_df[(kmer_df[freq_cols] > cons_thr).any(axis=1)]
    return out_kmer_df


def calc_kmer_pvalue(kmer: str, first_group, sec_group, matrix: pd.DataFrame):
    first_obs = (
        matrix.loc[kmer, first_group]
        .value_counts(sort=False)
        .reindex([0, 1], fill_value=0)
        .to_numpy()
    )
    sec_obs = (
        matrix.loc[kmer, sec_group]
        .value_counts(sort=False)
        .reindex([0, 1], fill_value=0)
        .to_numpy()
    )
    cont_table = np.concatenate([first_obs, sec_obs]).reshape(2, 2).T
    test_result = stats.fisher_exact(cont_table)
    return test_result.pvalue


def assign_kmer_group_ava(row: pd.Series, p_cols: list, avg_cols: list):
    max_group = row[avg_cols].idxmax()
    max_group = max_group.split("_")[0]
    rel_p_cols = [col for col in p_cols if max_group in col]
    if all(row[rel_p_cols] < 0.01):
        return f"{max_group}_enriched"
    else:
        return "not_significant"

def assign_kmer_group_ovr(row: pd.Series, goi:str, p_cols: list, err_cols: list):
   
    if all(row[p_cols] < 0.01) and all(row[err_cols] > 0):
        return f"{goi}_enriched"
    elif all(row[p_cols] < 0.01) and all(row[err_cols] < 0):
        return "control_enriched"
    elif any(row[p_cols] < 0.01):
        return "group_dependent"
    else:
        return "not_significant"

def check_exclusivity(row: pd.Series, avg_cols: list) -> str:
    if len([col for col in avg_cols if row[col] > 0]) == 1:
        return "exclusive"
    else:
        return "non-exclusive"


def build_similarity_graph_two_list(
    query_list, target_list, mismatch_treshold: int
) -> nx.Graph:

    tot = len(query_list)
    similarity_edges = set()
    i = 1
    t0 = datetime.now()
    for query_kmer in query_list:
        similarity_edges.add((query_kmer, query_kmer))
        for target_kmer in target_list:
            if (
                sum(c1 != c2 for c1, c2 in zip(query_kmer, target_kmer))
                <= mismatch_treshold
            ):
                similarity_edges.add((query_kmer, target_kmer))
        print(
            f"Done {i} out of {tot} kmers. Time per kmer: {datetime.now()-t0}",
            end="\r",
            flush=True,
        )
        t0 = datetime.now()
        i += 1

    silmilarity_graph = nx.Graph(similarity_edges)
    return silmilarity_graph


def calc_agg_freq(kmer_list, sample_list, source_df):
    pos_samples = set()
    for sample in sample_list:
        for kmer in kmer_list:
            if source_df.loc[kmer, sample] > 0:
                pos_samples.add(sample)
                break
    agg_freq = len(pos_samples) / len(sample_list)
    return agg_freq


# Deprecated
def map_kmers_find_mutations(kmer, ref_seq_str, pos_matrix, n=2, l=1000, find_wt=False):
    for i in range(n + 1):
        m = regex.findall(f"({kmer}){{s<={i}}}", ref_seq_str)
        if len(m) != 0:
            break
    if len(m) != 0:
        idx = [ref_seq_str.index(seq) + 1 for seq in m]
        mean_pos = pos_matrix.loc[kmer].mean()
        offset = list(abs(idx - mean_pos))

        if min(offset) <= l:
            best_align = offset.index(min(offset))
            best_match = m[best_align]
            start = idx[best_align]
            matchlist = [c1 == c2 for c1, c2 in zip(kmer, best_match)]
            n_match = sum(matchlist)
            alignment = {"seq": kmer, "start": start, "n_match": n_match}

            mutations = {"id": [], "ref_nt": [], "pos": [], "mut_nt": [], "kmer": []}
            for i, is_match in enumerate(matchlist):
                if find_wt == False:
                    if is_match == False:
                        mutations["ref_nt"].append(best_match[i])
                        mutations["pos"].append(start + i)
                        mutations["mut_nt"].append(kmer[i])
                        mutations["kmer"].append(kmer)
                        mutations["id"].append(f"{best_match[i]}{start+i}{kmer[i]}")

                if find_wt == True:
                    if is_match == True:
                        mutations["ref_nt"].append(best_match[i])
                        mutations["pos"].append(start + i)
                        mutations["mut_nt"].append(kmer[i])
                        mutations["kmer"].append(kmer)
                        mutations["id"].append(f"{best_match[i]}{start+i}{kmer[i]}")

        else:
            alignment = None
            mutations = {"id": [], "ref_nt": [], "pos": [], "mut_nt": [], "kmer": []}

    else:
        alignment = None
        mutations = {"id": [], "ref_nt": [], "pos": [], "mut_nt": [], "kmer": []}
    return alignment, mutations


def assemble_kmers(
    kmer_list: list, how: str = "pos", kmer_df: pd.DataFrame = None
) -> nx.DiGraph:
    edge_tuples = set()
    if how == "seq":
        for i in range(1, len(kmer_list)):
            for j in range(i):
                k1 = kmer_list[i]
                k2 = kmer_list[j]
                if k1[1:] == k2[:-1]:
                    edge_tuples.add((k1, k2))
                elif k2[1:] == k1[:-1]:
                    edge_tuples.add((k2, k1))
    elif how == "pos":
        if type(kmer_df) != pd.DataFrame:
            raise ValueError(
                "If assemblying kmers by position, you must provide a dataframe with kemrs as index and 'ref_start' column"
            )
        else:
            kmer_pos = [kmer_df.loc[kmer, "ref_start"] for kmer in kmer_list]
            for i in range(1, len(kmer_list)):
                for j in range(i):
                    k1 = kmer_list[i]
                    k2 = kmer_list[j]
                    if kmer_pos[i] == kmer_pos[j] - 1 and k1[1:] == k2[:-1]:
                        edge_tuples.add((k1, k2))
                    elif kmer_pos[j] == kmer_pos[i] - 1 and k2[1:] == k1[:-1]:
                        edge_tuples.add((k2, k1))

    kmer_graph = nx.from_edgelist(edge_tuples, create_using=nx.DiGraph)
    return kmer_graph


def plot_segments(segment_df, ref_seq, colormap=colormaps["coolwarm"]):
    plot_df = segment_df[["ref_start", "ref_end", "group"]].sort_values("ref_start")
    groups_values = ["ref"] + np.unique(plot_df["group"]).tolist()
    colors = [colormap(i) for i in np.linspace(0, 1, len(groups_values))]
    cdict = {groups_values[i]: colors[i] for i in range(len(groups_values))}
    ends_list = list(zip(plot_df["ref_start"], plot_df["ref_end"]))
    ys = []
    y = 1
    end_dict = {1: 0}
    while len(ends_list) > 1:
        ends = ends_list.pop(0)
        free_ends = {endy: endx for endy, endx in end_dict.items() if endx <= ends[0]}
        if len(free_ends.keys()) > 0:
            y = min(free_ends.keys())
        else:
            y = max(end_dict.keys()) + 1
        ys.append(y)
        end_dict[y] = ends[1]

    segments = [((1, 0), (len(ref_seq), 0))] + list(
        zip(zip(plot_df["ref_start"], ys), zip(plot_df["ref_end"], ys))
    )
    segment_colors = [cdict["ref"]] + [cdict[group] for group in plot_df["group"]]
    lines = collections.LineCollection(segments, colors=segment_colors, linewidths=2)
    fig, ax = pyplot.subplots()
    ax.set_xlim(0, len(ref_seq) + 1)
    ax.set_ylim(-1, max(ys) + 1)
    ax.add_collection(lines)
    legends = [patches.Patch(color=cdict[group], label=group) for group in cdict.keys()]
    pyplot.legend(handles=legends).set_loc("right")
    pyplot.show()


def _draw_genome_overlay_plot(
    segment_coords: list,
    segment_colors: list,
    ref_seq: str,
    title: str = None,
    colormap: colors.LinearSegmentedColormap = colormaps["coolwarm"],
    save_path: str = None,
    save_format: str = "svg"
):
    groups_values = ["ref"] + sorted(list(set(segment_colors)))
    colors = [colormap(i) for i in np.linspace(0, 1, len(groups_values))]
    cdict = {groups_values[i]: colors[i] for i in range(len(groups_values))}
    segment_colors = [cdict["ref"]] + [cdict[group] for group in segment_colors]
    segment_coords = [((1, 0), (len(ref_seq), 0))] + segment_coords
    lines = collections.LineCollection(
        segment_coords, colors=segment_colors, linewidths=2
    )
    fig, ax = pyplot.subplots()
    pyplot.tight_layout()
    ax.add_collection(lines)
    ax.autoscale()
    pyplot.tick_params(labelleft=False, left=False)
    ax.set_xlabel("Reference genome position")
    ax.set_ylabel("Genomes")
    ax.set_title(title)
    legends = [patches.Patch(color=cdict[group], label=group) for group in cdict.keys()]
    legends.reverse()
    pyplot.legend(handles=legends, bbox_to_anchor=(1.0, 0.5), loc="center left")
    if save_path != None:
        pyplot.savefig(save_path, format=save_format, dpi=600, bbox_inches="tight")
    pyplot.show()


def plot_kmers_by_genome(
    kmers_and_genomes:list,
    kmer_matrix:pd.DataFrame,
    group_sample_id:dict,
    ref_seq:str,
    title: str = None,
    colormap: colors.LinearSegmentedColormap = colormaps["coolwarm"],
    save_path: str = None,
    save_format: str = "svg"
):
    segment_coords = []
    segment_colors = []
    y = 1
    for kmers,samples in kmers_and_genomes:
        samples_to_plot = samples
        kmers_to_plot = kmers
        for sample_name in samples_to_plot:
            sample_group = [
                key for key, value in group_sample_id.items() if sample_name in value
            ][0]
            for kmer in kmers_to_plot:
                if kmer_matrix.loc[kmer, sample_name] > 0:
                    temp_segments = kmer_matrix.loc[kmer, "ref_pos"]
                    temp_segments = [((x1, y), (x2, y)) for x1, x2 in temp_segments]
                    segment_coords.extend(temp_segments)
                    segment_colors.extend([sample_group for _ in range(len(temp_segments))])
            y += 1
            
    _draw_genome_overlay_plot(segment_coords, segment_colors,ref_seq, title, colormap, save_path, save_format)

def assign_gene_from_interval(ref_pos: list, gene_dict: dict) -> str:
    genes = []
    for start, end in ref_pos:
        for gene, gene_ranges in gene_dict.items():
            for gene_range in gene_ranges:
                if (
                    len(
                        [
                            pos
                            for pos in range(start, end + 1)
                            if pos in list(range(gene_range[0], gene_range[1] + 1))
                        ]
                    )
                    > 0
                ):
                    genes.append(gene)

def assign_gene_from_position(ref_pos: int, gene_dict: dict) -> str:
    genes = []
    for gene, gene_ranges in gene_dict.items():
        for gene_range in gene_ranges:
            if gene_range[0] < ref_pos < gene_range[1]:
                genes.append(gene)
    return ", ".join(genes)


class AnalysisPipeline: 
    """
    AnalysisPipeline:
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
            all_kmer_matrix.index = all_kmer_matrix.index.map(lambda id: decode_kmer(id,k))

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
            common_kmer_matrix = filter_kmers(
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