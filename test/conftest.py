import pytest
import os
import shutil
from unittest import mock


@pytest.fixture
def correct_project_dir():
    cor_dir = "test/testproject"
    template_config = "templates/config.yaml"
    os.makedirs(cor_dir)
    shutil.copy2(template_config, cor_dir)
    yield cor_dir
    shutil.rmtree(f"{cor_dir}")


@pytest.fixture
def no_project_dir():
    fake_dir = "test/testproject"
    return fake_dir


@pytest.fixture
def empty_project_dir():
    empt_dir = "test/testproject"
    os.makedirs(empt_dir)
    yield empt_dir
    os.removedirs(empt_dir)


@pytest.fixture
def default_max_k():
    return 31


@pytest.fixture
def correct_k():
    return 21


@pytest.fixture
def float_k():
    return 31.0


@pytest.fixture
def small_k():
    return 3


@pytest.fixture
def even_k():
    return 30
