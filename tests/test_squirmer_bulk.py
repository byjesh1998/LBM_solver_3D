"""Single squirmer in a fully periodic box.

Checks
  * a neutral squirmer swims at U0 = 2 B1 / 3 (within the finite-box correction)
  * it swims straight: no lateral drift, no rotation
  * finite-Re trend: pusher faster, puller slower than neutral
  * near-field flow (1.25R < r < 3R, plane through the centre) agrees with the
    unbounded Blake solution to within 30 %. This is a regression check, not a
    precision test: R = 4 lattice units (staircase sphere), a 40^3 periodic box
    (image flows) and Re ~ 0.5 (Blake is a Stokes solution) give 15-25 %.
  * mass is conserved (with the global mass correction)
Figure: tests/output/squirmer_model.png
"""
import numpy as np
import pytest

from conftest import ROOT

INPUT = ROOT / "inputs" / "bulk_squirmer.in"
R, B1 = 4.0, 0.015
U0 = 2 * B1 / 3
CASES = [(-3, "pusher β = −3"), (0, "neutral β = 0"), (3, "puller β = +3")]
BOX, STEPS = 40, 2000


@pytest.fixture(scope="module")
def bulk_runs(lp, outdir):
    runs = {}
    for beta, label in CASES:
        prefix = f"bulk_beta{beta:+d}"
        lp.run(INPUT, dict(nx=BOX, ny=BOX, nz=BOX, beta=beta, steps=STEPS, radius=R, b1=B1, prefix=prefix),
               output_dir=outdir)
        runs[beta] = dict(
            label=label, beta=beta,
            traj=lp.read_traj(outdir / f"{prefix}_traj.dat"),
            state=lp.read_keyvalue(outdir / f"{prefix}_squirmer.dat"),
            summary=lp.read_keyvalue(outdir / f"{prefix}_summary.txt"),
            vtk=lp.read_vtk(outdir / f"{prefix}_final.vtk"),
        )
    return runs


def late(tr, key):
    return tr[key][tr["t"] > STEPS / 2]


def test_neutral_swimming_speed(bulk_runs):
    ratio = late(bulk_runs[0]["traj"], "Ux").mean() / U0
    assert 0.93 < ratio < 1.03, f"<Ux>/U0 = {ratio:.3f}"


def test_straight_swimming(bulk_runs):
    for beta, r in bulk_runs.items():
        tr = r["traj"]
        assert np.abs(late(tr, "Uy")).max() < 1e-3 * U0, f"beta={beta}: lateral drift"
        assert np.abs(late(tr, "Wz")).max() < 1e-6, f"beta={beta}: rotation"
        assert np.abs(tr["angle"]).max() < 1e-6


def test_finite_reynolds_speed_ordering(bulk_runs):
    u = {b: late(r["traj"], "Ux").mean() for b, r in bulk_runs.items()}
    assert u[-3] > u[0] > u[3], f"speeds {u}"


def test_mass_conservation(bulk_runs):
    for beta, r in bulk_runs.items():
        assert abs(r["summary"]["mean_density"] - 1.0) < 1e-4


def test_flow_field_matches_blake(lp, outdir, bulk_runs):
    cases = [dict(beta=b, label=r["label"], vtk=r["vtk"], state=r["state"]) for b, r in bulk_runs.items()]
    fig, errors = lp.plot_squirmer_model(cases, R, B1, path=outdir / "squirmer_model.png")
    for label, e in errors.items():
        assert e["near_field"] < 0.30, f"{label}: near-field error {e['near_field']:.3f}"
    assert (outdir / "squirmer_model.png").exists()
