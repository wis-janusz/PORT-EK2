import pytest
import portek
import pandas as pd

@pytest.fixture
def good_k():
    return 11

@pytest.fixture
def bad_k():
    return 27

@pytest.fixture
def project_dir_ovr():
    return "projects/test_ovr"

@pytest.fixture
def project_dir_ava():
    return "projects/test_ava"

@pytest.fixture
def wrong_project_dir():
    return "input"

@pytest.fixture
def no_project_dir():
    return "bleble"

@pytest.fixture
def good_c():
    return 0.5

@pytest.fixture
def good_min_rmse():
    return 0.5

def test_ava_get(project_dir_ava, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ava
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    test_pipeline.get_kmers(save_rare=True)
    assert type(test_pipeline.matrices["common"]) == pd.DataFrame
    assert type(test_pipeline.matrices["rare"]) == pd.DataFrame
    assert len(test_pipeline.matrices["rare"]) > len(test_pipeline.matrices["common"])


def test_ava_calc(project_dir_ava, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ava
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    test_pipeline.get_kmers(save_rare=True)
    test_pipeline.calc_kmer_stats("common")
    assert len(test_pipeline.enriched_groups) == 3
    assert "group" in test_pipeline.matrices["common"].columns


def test_ovr_get(project_dir_ovr, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ovr
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    test_pipeline.get_kmers(save_rare=True)
    assert type(test_pipeline.matrices["common"]) == pd.DataFrame
    assert type(test_pipeline.matrices["rare"]) == pd.DataFrame
    assert len(test_pipeline.matrices["rare"]) > len(test_pipeline.matrices["common"])


def test_ovr_calc(project_dir_ovr, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ovr
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    test_pipeline.get_kmers(save_rare=True)
    test_pipeline.calc_kmer_stats("common")
    assert len(test_pipeline.enriched_groups) == 2
    assert "group" in test_pipeline.matrices["common"].columns


def test_wrong_k(project_dir_ava, bad_k, good_c, good_min_rmse):
    proj_dir = project_dir_ava
    k = bad_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    with pytest.raises(FileNotFoundError):
        test_pipeline.get_kmers()
