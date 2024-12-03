import pytest
import portek
from unittest import mock

def test_KmerFinder_correct_dir(correct_project_dir, default_max_k):
    test_kmer_finder = portek.KmerFinder(correct_project_dir, default_max_k)
    assert test_kmer_finder.project_dir == "test/testproject/"
    assert test_kmer_finder.maxk == 31


def test_KmerFinder_no_dir(no_project_dir, default_max_k):
    with pytest.raises(NotADirectoryError):
        test_kmer_finder = portek.KmerFinder(no_project_dir, default_max_k)


def test_KmerFinder_empty_dir(empty_project_dir, default_max_k):
    with pytest.raises(FileNotFoundError):
        test_kmer_finder = portek.KmerFinder(empty_project_dir, default_max_k)


def test_KmerFinder_correct_k(correct_project_dir,correct_k):
    test_kmer_finder = portek.KmerFinder(correct_project_dir, correct_k)
    assert test_kmer_finder.maxk == 21


def test_KmerFinder_float_k(correct_project_dir,float_k):
    with pytest.raises(TypeError):
        test_kmer_finder = portek.KmerFinder(correct_project_dir, float_k)


def test_KmerFinder_small_k(correct_project_dir,small_k):
    with pytest.raises(TypeError):
        test_kmer_finder = portek.KmerFinder(correct_project_dir, small_k)


def test_KmerFinder_even_k(correct_project_dir,even_k):
    with pytest.raises(TypeError):
        test_kmer_finder = portek.KmerFinder(correct_project_dir, even_k)


correct_config_ava = """
sample_groups: [A,B,C] 
input_files: [A.fasta, B.fasta, C.fasta] 
header_format: [] 
mode: ava
goi: 
"""
@mock.patch("builtins.open", new_callable=mock.mock_open, read_data=correct_config_ava)
def test_KmerFinder_correct_config_ava(mock_config, correct_project_dir, default_max_k):
    test_kmer_finder = portek.KmerFinder(correct_project_dir, default_max_k)
    assert test_kmer_finder.project_dir == "test/testproject/"
    assert test_kmer_finder.maxk == 31
    assert test_kmer_finder.sample_groups == ["A", "B", "C"]