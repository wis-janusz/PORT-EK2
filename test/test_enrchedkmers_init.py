import pytest
import portek

@pytest.fixture
def good_k():
    return 11

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
    return 0.8

@pytest.fixture
def good_min_rmse():
    return 0.5

def test_genes(project_dir_ovr, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ovr
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    assert type(test_pipeline.ref_genes == dict) 
    assert type(test_pipeline.ref_genes["5LTR"][0] == list)
    assert test_pipeline.ref_genes["5LTR"][0][0] == 1


def test_ava(project_dir_ava, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ava
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    assert test_pipeline.sample_groups == ["MD","N","O"]
    assert test_pipeline.mode == "ava"
    assert test_pipeline.goi == None
    assert test_pipeline.control_groups == None


def test_ovr(project_dir_ovr, good_k, good_c, good_min_rmse):
    proj_dir = project_dir_ovr
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)
    assert test_pipeline.sample_groups == ["MD","N","O"]
    assert test_pipeline.mode == "ovr"
    assert test_pipeline.goi == "MD"
    assert test_pipeline.control_groups == ["N","O"]


def test_no_config(wrong_project_dir, good_k, good_c, good_min_rmse):
    proj_dir = wrong_project_dir
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    with pytest.raises(FileNotFoundError):
        test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)

def test_no_dir(no_project_dir, good_k, good_c, good_min_rmse):
    proj_dir = no_project_dir
    k = good_k
    c = good_c
    min_rmse = good_min_rmse
    with pytest.raises(NotADirectoryError):
        test_pipeline = portek.EnrichedKmersPipeline(proj_dir, k, c, min_rmse)