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