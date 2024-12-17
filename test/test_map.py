import pytest
import portek
from unittest import mock

import pandas as pd
import numpy as np

def test_check_bowtie2_exists(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    assert mapper._check_bowtie2_path() == True


def test_check_bowtie2_notexists(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    mapper.bowtie2_path = "test/testproject"
    assert mapper._check_bowtie2_path() == False


def test_check_index_exists(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    assert mapper._check_index_built() == True


def test_check_index_notexists(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    mapper.ref_seq_name = "foo"
    assert mapper._check_index_built() == False


def test_bowtie_build_index_correct(
    correct_project_dir, correct_k, correct_bowtie_build_cmd
):
    mock_result = mock.Mock()
    mock_result.returncode = 0
    with mock.patch("subprocess.run", return_value=mock_result) as mock_subprocess:
        mapper = portek.MappingPipeline(correct_project_dir, correct_k)
        result = mapper._bowtie_build_index(verbose=False)
        expected_string = correct_bowtie_build_cmd
        call_args = mock_subprocess.call_args
        assert expected_string == " ".join(*call_args[0])


def test_bowtie_map_correct(correct_project_dir, correct_k, correct_bowtie_map_cmd):
    mock_result = mock.Mock()
    mock_result.returncode = 0
    with mock.patch("subprocess.run", return_value=mock_result) as mock_subprocess:
        mapper = portek.MappingPipeline(correct_project_dir, correct_k)
        result = mapper._bowtie_map(verbose=False)
        expected_string = correct_bowtie_map_cmd
        call_args = mock_subprocess.call_args
        assert expected_string == " ".join(*call_args[0])


def test_read_sam_correct(correct_project_dir, correct_k, test_mapping_groups):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    mapping_df = mapper._read_sam_to_df()
    assert len(mapping_df) == 65
    assert mapping_df.columns.equals(
        pd.Index(["kmer", "flag", "ref_pos", "CIGAR", "n_mismatch", "group"])
    )
    assert mapping_df["group"].to_list() == test_mapping_groups


def test_parse_CIGAR(correct_project_dir, correct_k, CIGARS_to_parse, expected_CIGARS):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    parsed_CIGARS = []
    for cigar in CIGARS_to_parse:
        parsed_CIGARS.append(mapper._parse_CIGAR(cigar))
    assert parsed_CIGARS == expected_CIGARS


def test_detect_CIGAR_unmapped(correct_project_dir, correct_k, expected_CIGARS):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    clippings = []
    for cigar in expected_CIGARS:
        clippings.append(mapper._detect_unmapped_CIGAR(cigar))
    assert clippings == [False, False, False, True, True, True]


def test_align_seqs(correct_project_dir, correct_k, test_ref_seq, test_mappings):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    expected_aln_lens = [7, 7, 7, 10, 13]
    expected_aln_q = [
        ["C", "T", "G", "T", "C", "G", "C"],
        ["C", "T", "C", "T", "C", "G", "C"],
        ["C", "C", "C", "A", "A", "G", "C"],
        ["C", "C", "C", "-", "-", "-", "G", "T", "C", "G"],
        ["C", "C", "C", "A", "A", "G", "C", "-", "-", "-", "C", "C", "C"],
    ]
    expected_aln_t = [
        ["C", "T", "G", "T", "C", "G", "C"],
        ["C", "T", "G", "T", "C", "G", "C"],
        ["C", "C", "C", "-", "-", "G", "C"],
        ["C", "C", "C", "G", "C", "T", "G", "T", "C", "G"],
        ["C", "C", "C", "-", "-", "G", "C", "T", "G", "T", "C", "G", "C"],
    ]
    expected_aln_pos = [
        [5, 6, 7, 8, 9, 10, 11],
        [5, 6, 7, 8, 9, 10, 11],
        [1, 2, 3, 3, 3, 4, 5],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    ]

    for i in range(5):
        aligned_mapping = mapper._align_seqs(
            test_ref_seq,
            test_mappings["kmer"][i],
            test_mappings["pos"][i],
            test_mappings["CIGAR"][i],
        )
        assert (
            len(aligned_mapping[0])
            == len(aligned_mapping[1])
            == len(aligned_mapping[2])
            == expected_aln_lens[i]
        )
        assert aligned_mapping[0] == expected_aln_q[i]
        assert aligned_mapping[1] == expected_aln_t[i]        
        assert aligned_mapping[2] == expected_aln_pos[i]


def test_join_indels(
    correct_project_dir, correct_k, expected_mutations, expected_mutations_joined
):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    for i in range(5):
        mutations = mapper._join_indels(expected_mutations[i])
        assert mutations == expected_mutations_joined[i]


def test_find_variants(
    correct_project_dir,
    correct_k,
    test_ref_seq,
    test_mappings,
    expected_mutations_joined,
):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    for i in range(5):
        mutations = mapper._find_variants(
            test_ref_seq,
            test_mappings["kmer"][i],
            test_mappings["pos"][i],
            test_mappings["CIGAR"][i],
        )
        assert mutations == expected_mutations_joined[i]


def test_mutations_tuples_to_text(
    correct_project_dir, correct_k, expected_mutations_joined, expected_mutations_text
):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    for i in range(5):
        mutations = mapper._mutation_tuples_to_text(expected_mutations_joined[i])
        assert mutations == expected_mutations_text[i]
