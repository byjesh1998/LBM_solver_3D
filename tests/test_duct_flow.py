"""Fluid solver without a squirmer: Poiseuille flow in a square duct and a slit.

Checks
  * relative L2 error vs the exact duct solution (series) and plane Poiseuille flow
  * second-order grid convergence at fixed Reynolds and Mach numbers
  * global momentum balance: force on the walls = total body force
  * mass conservation
Figure: tests/output/duct_validation.png
"""
import numpy as np
import pytest

from conftest import ROOT

INPUT = ROOT / "inputs" / "duct.in"
TAU, NU = 0.8, (0.8 - 0.5) / 3
WIDTHS = [8, 16, 32]


def run_duct(lp, outdir, H, **extra):
    fx = 5e-5 * (32 / H) ** 2              # keeps u_max (hence Re and Ma) fixed
    ov = dict(nx=4, ny=H + 2, nz=H + 2, tau=TAU, fx=fx, prefix=f"duct_H{H}", **extra)
    lp.run(INPUT, ov, output_dir=outdir)
    summary = lp.read_keyvalue(outdir / f"duct_H{H}_summary.txt")
    return fx, summary


@pytest.fixture(scope="module")
def duct_runs(lp, outdir):
    return {H: run_duct(lp, outdir, H) for H in WIDTHS}


def test_duct_matches_exact_solution(duct_runs):
    _, s = duct_runs[32]
    assert s["converged"] == 1
    assert s["l2_error"] < 5e-4, f"L2 error {s['l2_error']:.2e} too large"


def test_second_order_convergence(duct_runs):
    err = np.array([duct_runs[H][1]["l2_error"] for H in WIDTHS])
    order = np.log(err[:-1] / err[1:]) / np.log(2)
    assert np.all(order > 1.8), f"observed orders {order}"


def test_momentum_balance_and_mass(duct_runs):
    for H, (_, s) in duct_runs.items():
        fw, fb = s["wall_force"][0], s["body_force_total"][0]
        assert abs(fw - fb) / fb < 1e-6, f"H={H}: wall force {fw} vs body force {fb}"
        assert abs(s["mean_density"] - 1.0) < 1e-10


def test_plane_channel(lp, outdir):
    """Walls in y only: TRT with Lambda = 3/16 is exact for plane Poiseuille flow."""
    lp.run(INPUT, dict(nx=4, ny=34, nz=1, walls_z=0, tau=0.8, fx=1e-5, prefix="plane"), output_dir=outdir)
    s = lp.read_keyvalue(outdir / "plane_summary.txt")
    assert s["l2_error"] < 1e-6


def test_duct_figure(lp, outdir, duct_runs):
    fx, _ = duct_runs[32]
    vtk = lp.read_vtk(outdir / "duct_H32_final.vtk")
    errors = [duct_runs[H][1]["l2_error"] for H in WIDTHS]
    fig, err = lp.plot_duct_validation(vtk, fx, NU, convergence=(WIDTHS, errors),
                                       path=outdir / "duct_validation.png")
    assert abs(err - duct_runs[32][1]["l2_error"]) < 1e-8   # python and C++ agree
    assert (outdir / "duct_validation.png").exists()
