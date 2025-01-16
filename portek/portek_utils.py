import numpy as np
import pandas as pd
import networkx as nx
import regex
from matplotlib import pyplot, collections, colormaps, patches, colors
from datetime import datetime
from scipy import stats
from scipy.spatial import distance
from scipy.cluster import hierarchy
from Bio import SeqIO, SeqRecord, Seq


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


def calc_kmer_pvalue(kmer: str, first_group, sec_group, matrix: pd.DataFrame, freq_cols, err_cols):
    # if all(matrix.loc[kmer, freq_cols] > 0.9) and all(abs(matrix.loc[kmer, err_cols]) < 0.1):
    #     return 1
    # else:
    first_obs = matrix.loc[kmer, first_group]
    sec_obs = matrix.loc[kmer, sec_group]
    test_result = stats.mannwhitneyu(first_obs, sec_obs)
    return test_result.pvalue


def assign_kmer_group_ava(
    row: pd.Series, p_cols: list, avg_cols: list, freq_cols: list, err_cols: list
):
    max_group = row[avg_cols].idxmax()
    max_group = max_group.split("_")[0]
    rel_p_cols = [col for col in p_cols if max_group in col]
    if all(row[freq_cols] > 0.9) and all((abs(row[err_cols]) < 0.1)):
        return "conserved"
    elif all(row[rel_p_cols] < 0.01):
        return f"{max_group}_enriched"
    else:
        return "not_significant"


def assign_kmer_group_ovr(
    row: pd.Series, goi: str, p_cols: list, err_cols: list, freq_cols: list
):
    if all(row[freq_cols] > 0.9) and all((abs(row[err_cols]) < 0.1)):
        return "conserved"
    elif all(row[p_cols] < 0.01) and all(row[err_cols] > 0):
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

#Deprecated
def build_similarity_graph_two_list(
    name, query_list, target_list, mismatch_treshold: int
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
    return name, silmilarity_graph

#Deprecated
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


def assemble_kmers_debruijn(kmers: list) -> nx.MultiDiGraph:
    de_bruijn = nx.MultiDiGraph()
    for kmer in kmers:
        kmerL = kmer[:-1]
        kmerR = kmer[1:]
        de_bruijn.add_node(kmerL)
        de_bruijn.add_node(kmerR)
        de_bruijn.add_edge(kmerL, kmerR, sequence=kmer)
    return de_bruijn


def assemble_contig(contig_graph: nx.MultiDiGraph) -> str:
    kmers = []
    if contig_graph.number_of_edges() == 1:
        path = nx.eulerian_path(contig_graph, keys=True)
        for edge in path:
            sequence = "".join([edge[0], edge[1][-1:]])
            kmers.append(nx.get_edge_attributes(contig_graph, "sequence")[edge])
    else:
        path = nx.eulerian_path(contig_graph, keys=True)
        first = True
        for edge in path:
            if first == True:
                sequence = "".join([edge[0][:1], edge[1]])
                first = False
            else:
                sequence = "".join([sequence, edge[1][-1:]])
            kmers.append(nx.get_edge_attributes(contig_graph, "sequence")[edge])
    return sequence, kmers


def cluster_kmers(matrix: pd.DataFrame, discard_singles:bool = False, max_d:float = 0):
    distances = distance.pdist(matrix)
    linkage = hierarchy.linkage(distances)
    clustering = hierarchy.fcluster(linkage, max_d, criterion="distance")
    clustering = pd.Series(clustering, index=matrix.index, name="freq_cluster")
    freq_clusters = {cluster_no: [] for cluster_no in clustering.unique()}
    for cluster_no in freq_clusters.keys():
        freq_clusters[cluster_no] = clustering.loc[clustering == cluster_no].index.to_list()
    cluster_graphs = {}
    for cluster, kmers in freq_clusters.items():
        cluster_graphs[cluster] = assemble_kmers_debruijn(kmers)
    contigs_kmer_dict = {}
    kmer_contig_dict = {}
    for cluster, graph in cluster_graphs.items():
        components = [
            graph.subgraph(component).copy()
            for component in nx.weakly_connected_components(graph)
        ]
        for component in components:
            if nx.is_directed_acyclic_graph(component):
                if nx.has_eulerian_path(component):
                    if discard_singles == True:
                        if component.number_of_edges() == 1:
                            continue
                    sequence, kmer_list = assemble_contig(component)
                    contigs_kmer_dict[sequence] = kmer_list
                    for kmer in kmer_list:
                        kmer_contig_dict[kmer] = sequence
                else:
                    print(f"Ambiguous assembly in cluster {cluster}")
            else:
                print(f"Ambiguous assembly in cluster {cluster}")
    return kmer_contig_dict, contigs_kmer_dict

#deprecated
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

#deprecated
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

#deprecated
def _draw_genome_overlay_plot(
    segment_coords: list,
    segment_colors: list,
    ref_seq: str,
    title: str = None,
    colormap: colors.LinearSegmentedColormap = colormaps["coolwarm"],
    save_path: str = None,
    save_format: str = "svg",
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

#deprecated
def plot_kmers_by_genome(
    kmers_and_genomes: list,
    kmer_matrix: pd.DataFrame,
    group_sample_id: dict,
    ref_seq: str,
    title: str = None,
    colormap: colors.LinearSegmentedColormap = colormaps["coolwarm"],
    save_path: str = None,
    save_format: str = "svg",
):
    segment_coords = []
    segment_colors = []
    y = 1
    for kmers, samples in kmers_and_genomes:
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
                    segment_colors.extend(
                        [sample_group for _ in range(len(temp_segments))]
                    )
            y += 1

    _draw_genome_overlay_plot(
        segment_coords, segment_colors, ref_seq, title, colormap, save_path, save_format
    )

#deprecated
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

#deprecated
def assign_gene_from_position(ref_pos: int, gene_dict: dict) -> str:
    genes = []
    for gene, gene_ranges in gene_dict.items():
        for gene_range in gene_ranges:
            if gene_range[0] < ref_pos < gene_range[1]:
                genes.append(gene)
    return ", ".join(genes)


def make_ordinal(n):
    """
    Convert an integer into its ordinal representation::

        make_ordinal(0)   => '0th'
        make_ordinal(3)   => '3rd'
        make_ordinal(122) => '122nd'
        make_ordinal(213) => '213th'
    """
    n = int(n)
    if 11 <= (n % 100) <= 13:
        suffix = "th"
    else:
        suffix = ["th", "st", "nd", "rd", "th"][min(n % 10, 4)]
    return str(n) + suffix


def save_kmers_fasta(kmers:list, ids:list | None, name:str, directory:str, k:int) -> None:
    out_fasta_list = []
    for i,kmer in enumerate(kmers):
        if ids == None:
            out_fasta_list.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(kmer),
                    id=f"{kmer}",
                    description="",
                )
            )
        else:
            out_fasta_list.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(kmer),
                    id=ids[i],
                    description="",
                )
            )           
    SeqIO.write(
        out_fasta_list,
        f"{directory}/temp/{name}_{k}mers.fasta",
        format="fasta",
    )