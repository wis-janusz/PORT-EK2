import pytest
import os
import shutil
from unittest import mock


@pytest.fixture
def correct_project_dir():
    cor_dir = "test/testproject/"
    return cor_dir


@pytest.fixture
def no_project_dir():
    fake_dir = "test/faketestproject/"
    return fake_dir


@pytest.fixture
def empty_project_dir():
    empt_dir = "test/emptytestproject/"
    os.makedirs(empt_dir)
    yield empt_dir
    os.removedirs(empt_dir)


@pytest.fixture
def default_max_k():
    return 31


@pytest.fixture
def correct_k():
    return 15


@pytest.fixture
def float_k():
    return 31.0


@pytest.fixture
def small_k():
    return 3


@pytest.fixture
def even_k():
    return 30


@pytest.fixture
def correct_bowtie_build_cmd():
    return "test/mock_bowtie2//bowtie2-build -f test/testproject//input/ref_seq.fasta test/testproject//temp/ref_index/ref_seq"


@pytest.fixture
def correct_bowtie_map_cmd():
    return "test/mock_bowtie2//bowtie2 -a --norc --no-hd -L 8 -x test/testproject//temp/ref_index/ref_seq -f test/testproject//temp/enriched_15mers.fasta -S test/testproject//temp/enriched_15mers.sam"


@pytest.fixture
def CIGARS_to_parse():
    cigars = ["15M", "8M3I4M", "8M1D1M2I4M", "10M5H", "2S10M3S"]
    return cigars


@pytest.fixture
def expected_CIGARS():
    cigars = [
        {"M": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]},
        {"M": [0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14], "I": [8, 9, 10]},
        {"M": [0, 1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 14, 15], "D": [8], "I": [10, 11]},
        {"M": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], "H": [10, 11, 12, 13, 14]},
        {"S": [0, 1, 12, 13, 14], "M": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]},
    ]
    return cigars