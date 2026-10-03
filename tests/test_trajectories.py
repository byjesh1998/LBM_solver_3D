"""Optional (slow) tests: squirmer trajectories in a slit channel.

Run with  python -m pytest tests -v --runslow   (~15 min on one core)

Checks
  * neutral squirmer (beta = 0): oscillates between both walls
  * puller (beta = +1): the lateral oscillation is damped
Figure: tests/output/trajectories.png
"""
import numpy as np
import pytest

from conftest import ROOT

NY, R, B1 = 34, 4.0, 0.015
U0 = 2 * B1 / 3


@pytest.fixture(scope="module")
def channel_runs(lp, outdir):
    runs = {}
    for name, steps in [("neutral", 20000), ("puller", 42000)]:
        prefix = f"traj_{name}"
        lp.run(ROOT / "inputs" / f"channel_{name}.in", dict(steps=steps, prefix=prefix), output_dir=outdir)
        runs[name] = lp.read_traj(outdir / f"{prefix}_traj.dat")
    return runs


@pytest.mark.slow
def test_neutral_oscillates(lp, channel_runs):
    T, yn, idx = lp.oscillation_extrema(channel_runs["neutral"], NY, R, U0, prominence=0.3)
    assert yn.max() > 0.7 and yn.min() < -0.7, "neutral squirmer should visit both walls"
    assert np.any(np.diff(np.sign(yn[T > 1.5])) != 0), "should cross the centreline"


@pytest.mark.slow
def test_puller_oscillation_is_damped(lp, channel_runs):
    T, yn, idx = lp.oscillation_extrema(channel_runs["puller"], NY, R, U0, t_skip=1.5, prominence=0.05)
    assert len(idx) >= 2, "need two extrema after the first wall contact"
    amp = np.abs(yn[idx])
    assert amp[1] < amp[0] - 0.05, f"amplitudes {amp} should decrease"


@pytest.mark.slow
def test_trajectory_figure(lp, outdir, channel_runs):
    lp.plot_trajectories([(channel_runs["neutral"], "neutral β = 0", "k"),
                          (channel_runs["puller"], "puller β = +1", "C0")],
                         NY, R, U0, path=outdir / "trajectories.png")
    assert (outdir / "trajectories.png").exists()
