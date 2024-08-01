import numpy as np
import pandas as pd
import networkx as nx
import regex

from datetime import datetime
from scipy import stats


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


# Future upgrade: if results found for i, but more than l away and thus discarded, go back to i+1.
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


def map_kmers(
    enriched_kmers_df: pd.DataFrame,
    pos_matrix: pd.DataFrame,
    sample_list: list,
    ref_seq: str,
    m_map=2,
    l_map=1000,
) -> pd.DataFrame:
    
    non_zero_pos_matrix = pos_matrix[pos_matrix > 0]
    kmer_mapping = pd.DataFrame(0, index = enriched_kmers_df.index, columns = ["ref_pos","ref_n_match"])

    for kmer in enriched_kmers_df.index:
        if enriched_kmers_df.loc[kmer,sample_list].max(axis=1) <= 1:
            matches = []
            for i in range(m_map + 1):
                matches.append(regex.findall(f"({kmer}){{s<={i}}}", ref_seq))

            for matchlist in matches:
                if len(matchlist) != 0:
                    ref_pos = [ref_seq.index(seq) + 1 for seq in matchlist]
                    mean_pos = non_zero_pos_matrix.loc[kmer].mean()
                    offset = list(abs(ref_pos - mean_pos))

                    if min(offset) <= l_map:
                        best_align = offset.index(min(offset))
                        best_match = matchlist[best_align]
                        start = ref_pos[best_align]
                        nuc_match = [c1 == c2 for c1, c2 in zip(kmer, best_match)]
                        n_match = sum(nuc_match)
                        kmer_mapping.loc[kmer,"ref_pos"] = start
                        kmer_mapping.loc[kmer,"ref_n_match"] = n_match
                        break

    return kmer_mapping


def calc_nuc_distr(enriched_kmers_df: pd.DataFrame, ref_seq:str, freq_cols:list):
    distr_cols = []
    for col in freq_cols:
        group = col.replace("_freq", "")
        for nuc in ["A","C","G","T"]:
            distr_cols.append[f"{nuc}_{group}"]

    nuc_distr = pd.DataFrame(0, index=pd.RangeIndex(1,len(ref_seq)+1), columns=distr_cols)

    for kmer in enriched_kmers_df.index:
        if enriched_kmers_df.loc[kmer,"ref_pos"] != 0:
            for i, nuc in enumerate(kmer):
                nuc_distr
