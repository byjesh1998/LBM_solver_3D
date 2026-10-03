"""Shared pytest configuration.

Fast tests (default):   python -m pytest tests -v            (~1.5 min)
Including slow tests:   python -m pytest tests -v --runslow  (~15 min more)

Figures are written to tests/output/.
"""
import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "postprocessing"))
import lbm_post  # noqa: E402

OUTPUT = ROOT / "tests" / "output"


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False,
                     help="also run the long squirmer-trajectory tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long trajectory runs (enable with --runslow)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="slow test: use --runslow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope="session")
def lp():
    """The post-processing module, with the solver compiled."""
    lbm_post.build()
    return lbm_post


@pytest.fixture(scope="session")
def outdir():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    return OUTPUT
