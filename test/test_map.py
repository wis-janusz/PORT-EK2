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


def test_read_sam_correct(correct_project_dir, correct_k):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    mapping_df = mapper._read_sam_to_df()
    assert mapping_df.to_numpy().shape == (64, 3)
    assert len(mapping_df.dtypes.unique()) == 2


def test_parse_CIGAR(correct_project_dir, correct_k, CIGARS_to_parse, expected_CIGARS):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    parsed_CIGARS = []
    for cigar in CIGARS_to_parse:
        parsed_CIGARS.append(mapper._parse_CIGAR(cigar))
    assert parsed_CIGARS == expected_CIGARS


def test_detect_CIGAR_clip(correct_project_dir, correct_k, CIGARS_to_parse):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    clippings = []
    for cigar in CIGARS_to_parse:
        clippings.append(mapper._detect_clip_CIGAR(cigar))
    assert clippings == [False, False, False, True, True]


def test_align_seqs(correct_project_dir, correct_k, test_ref_seq, test_mappings):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    expected_aln_lens = [7, 7, 7, 10, 13]
    expected_aln_pos = [
        [6, 7, 8, 9, 10, 11, 12],
        [6, 7, 8, 9, 10, 11, 12],
        [2, 3, 4, 4, 4, 5, 6],
        [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        [2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12],
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
        assert (
            aligned_mapping[2]
            == expected_aln_pos[i]
        )



def test_find_variants(correct_project_dir, correct_k, test_ref_seq, test_mappings):
    mapper = portek.MappingPipeline(correct_project_dir, correct_k)
    expected_mutations = [
        [],
        ["8G>C"],
        ["4insA", "4insA"],
        ["5del", "6del", "7del"],
        ["4insA", "4insA", "7del", "8del", "9del", "11G>C"],
    ]
    for i in range(5):
        mutations = mapper._find_variants(
            test_ref_seq,
            test_mappings["kmer"][i],
            test_mappings["pos"][i],
            test_mappings["CIGAR"][i],
            test_mappings["n_mismatch"][i],
        )
        assert mutations == expected_mutations[i]
