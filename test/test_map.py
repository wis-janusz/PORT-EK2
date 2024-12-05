import pytest
import portek
from unittest import mock


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
    mapper.ref_seq = "foo"
    assert mapper._check_index_built() == False 


def test_bowtie_build_index_correct(correct_project_dir, correct_k, correct_bowtie_build_cmd):
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

def test_read_sam_correct(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    mapping_df = mapper._read_sam_to_df()
    assert mapping_df.to_numpy().shape == (64,3)
    assert len(mapping_df.dtypes.unique()) == 2

def test_parse_CIGAR(correct_project_dir, correct_k,CIGARS_to_parse, expected_CIGARS):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    parsed_CIGARS = []
    for cigar in CIGARS_to_parse:
        parsed_CIGARS.append(mapper._parse_CIGAR(cigar))
    assert parsed_CIGARS == expected_CIGARS

def test_analyze_mapping(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    test_df = mapper.analyze_mapping()
    assert test_df.to_numpy().shape == (63,4)
    assert len(test_df.dtypes.unique()) == 2
