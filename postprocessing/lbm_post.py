"""
lbm_post -- post-processing helpers for lbm-squirmer.

Running           : build(), run()
Reading output    : read_vtk(), read_traj(), read_keyvalue()
Exact solutions   : duct_velocity(), plane_velocity(), blake_velocity()
Analysis          : duct_cross_section_error(), squirmer_plane(), oscillation_extrema()
Plotting          : plot_duct_validation(), plot_squirmer_model(), plot_trajectories()

All quantities are in lattice units (dx = dt = 1, rho0 = 1).
"""


from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "bin" / "lbm_squirmer"


# ----------------------------------------------------------------------------
# Running the solver
# ----------------------------------------------------------------------------
def build(quiet: bool = True) -> Path:
    """Compile the solver with `make` (no-op if it is up to date)."""
    res = subprocess.run(["make", "-C", str(ROOT)], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError("build failed:\n" + res.stdout + res.stderr)
    if not quiet:
        print(res.stdout.strip() or "up to date")
    return BINARY


def run(input_file: str | os.PathLike, overrides: dict | None = None,
        output_dir: str | os.PathLike | None = None, echo: bool = False) -> str:
    """Run the solver on an input file with optional key=value overrides.

    Returns the solver's standard output. Raises RuntimeError on failure.
    """
    args = [str(BINARY), str(input_file)]
    overrides = dict(overrides or {})
    if output_dir is not None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        overrides["output_dir"] = str(output_dir)
    args += [f"{k}={v}" for k, v in overrides.items()]
    res = subprocess.run(args, capture_output=True, text=True, cwd=ROOT)
    if res.returncode != 0:
        raise RuntimeError(f"solver failed ({' '.join(args)}):\n{res.stdout}{res.stderr}")
    if echo:
        print(res.stdout)
    return res.stdout


# ----------------------------------------------------------------------------
# Reading output
# ----------------------------------------------------------------------------
_VTK_TYPES = {"unsigned_char": ">u1", "char": ">i1", "int": ">i4", "float": ">f4", "double": ">f8"}


def read_vtk(path) -> dict:
    """Read a legacy binary STRUCTURED_POINTS file written by the solver.

    Returns a dict with 'dims' (nx, ny, nz) and every field as an array
    indexed [z, y, x] (scalars) or [z, y, x, component] (vectors).
    """
    data = Path(path).read_bytes()
    pos = 0

    def line():
        nonlocal pos
        end = data.index(b"\n", pos)
        s = data[pos:end].decode().strip()
        pos = end + 1
        return s

    out, npts, dims = {}, None, None
    while pos < len(data):
        s = line()
        if not s:
            continue
        tok = s.split()
        if tok[0] == "DIMENSIONS":
            dims = tuple(int(v) for v in tok[1:4])
        elif tok[0] == "POINT_DATA":
            npts = int(tok[1])
        elif tok[0] in ("SCALARS", "VECTORS"):
            name, typ = tok[1], tok[2]
            ncomp = 3 if tok[0] == "VECTORS" else (int(tok[3]) if len(tok) > 3 else 1)
            if tok[0] == "SCALARS":
                line()  # LOOKUP_TABLE
            dt = np.dtype(_VTK_TYPES[typ])
            nbytes = npts * ncomp * dt.itemsize
            arr = np.frombuffer(data[pos:pos + nbytes], dtype=dt).astype(dt.newbyteorder("="))
            pos += nbytes
            nx, ny, nz = dims
            out[name] = arr.reshape(nz, ny, nx, ncomp) if ncomp > 1 else arr.reshape(nz, ny, nx)
    out["dims"] = dims
    return out


TRAJ_COLUMNS = ["t", "X", "Y", "Z", "ex", "ey", "ez", "Ux", "Uy", "Uz",
                "Wx", "Wy", "Wz", "Fwall_y", "Fwall_z"]


def read_traj(path) -> dict:
    """Read <prefix>_traj.dat into a dict of 1D arrays (see TRAJ_COLUMNS).

    Also adds 'angle' (orientation in the x-y plane, degrees).
    """
    d = np.atleast_2d(np.loadtxt(path))
    out = {name: d[:, k] for k, name in enumerate(TRAJ_COLUMNS[: d.shape[1]])}
    out["angle"] = np.degrees(np.arctan2(out["ey"], out["ex"]))
    return out


def read_keyvalue(path) -> dict:
    """Read 'key value [value ...]' text files (_summary.txt, _squirmer.dat)."""
    out = {}
    for s in Path(path).read_text().splitlines():
        tok = s.split()
        if not tok:
            continue
        vals = [float(v) for v in tok[1:]]
        out[tok[0]] = vals[0] if len(vals) == 1 else np.array(vals)
    return out


# ----------------------------------------------------------------------------
# Exact solutions
# ----------------------------------------------------------------------------
def duct_velocity(y, z, a, b, G, mu, nterms=200):
    """Rectangular duct -a<y<a, -b<z<b driven by force density G (series solution)."""
    y, z = np.asarray(y, float), np.asarray(z, float)
    s = np.zeros(np.broadcast(y, z).shape)
    for n in range(1, 2 * nterms, 2):
        k = n * np.pi / (2 * a)
        A, B = k * np.abs(z), k * b
        ratio = np.exp(A - B) * (1 + np.exp(-2 * A)) / (1 + np.exp(-2 * B))
        s = s + (-1) ** ((n - 1) // 2) * (1 - ratio) * np.cos(k * y) / n ** 3
    return 16 * a * a * G / (mu * np.pi ** 3) * s


def plane_velocity(y, a, G, mu):
    """Plane Poiseuille flow -a<y<a."""
    return G / (2 * mu) * (a * a - np.asarray(y, float) ** 2)


def blake_velocity(x, y, z, R, B1, beta, e=(1.0, 0.0, 0.0)):
    """Lab-frame Stokes flow around a squirmer centred at the origin (Blake 1971).

    Returns an array of shape (3, *x.shape); NaN inside the sphere.
    """
    x, y, z = np.broadcast_arrays(*(np.asarray(v, float) for v in (x, y, z)))
    e = np.asarray(e, float)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    r = np.maximum(r, 1e-12)
    rh = np.stack([x / r, y / r, z / r])
    ct = np.tensordot(e, rh, axes=1)
    E = e.reshape((3,) + (1,) * x.ndim)
    B2 = beta * B1
    u = (B1 * (R / r) ** 3 * (ct * rh - E / 3)
         + B2 * ((R / r) ** 4 - (R / r) ** 2) * (1.5 * ct ** 2 - 0.5) * rh
         + B2 * (R / r) ** 4 * ct * (ct * rh - E))
    u[:, r < R] = np.nan
    return u


# ----------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------
def duct_cross_section_error(vtk: dict, G: float, nu: float, walls_y=True, walls_z=True):
    """LBM u_x on the cross-section x=0, the exact solution and the relative L2 error.

    Returns (Y, Z, u_lbm, u_exact, err) with Y, Z measured from the channel centre.
    """
    nx, ny, nz = vtk["dims"]
    ux = vtk["velocity"][..., 0][:, :, 0]
    fl = vtk["flag"][:, :, 0]
    yc, zc = 0.5 * (ny - 1), 0.5 * (nz - 1)
    a, b = 0.5 * (ny - 2), 0.5 * (nz - 2)
    Y, Z = np.meshgrid(np.arange(ny) - yc, np.arange(nz) - zc)
    if walls_y and walls_z:
        ue = duct_velocity(Y, Z, a, b, G, nu)
    elif walls_y:
        ue = plane_velocity(Y, a, G, nu)
    else:
        ue = plane_velocity(Z, b, G, nu)
    fluid = fl == 0
    u = np.where(fluid, ux, np.nan)
    ue = np.where(fluid, ue, np.nan)
    err = np.sqrt(np.nansum((u - ue) ** 2) / np.nansum(ue ** 2))
    return Y, Z, u, ue, err


def squirmer_plane(vtk: dict, state: dict, periodic=(True, True, True)):
    """Velocity in the x-y lattice plane nearest to the squirmer centre, in
    coordinates relative to the centre (minimum image in periodic directions).

    Returns (px, py, dz, ux, uy) with px sorted ascending; NaN on solid nodes.
    """
    nx, ny, nz = vtk["dims"]
    X, Y, Z = state["X"]
    zp = int(np.floor(Z % nz)) if periodic[2] else int(np.floor(Z))
    dz = zp - (Z % nz if periodic[2] else Z)
    dx = np.arange(nx) - (X % nx)
    if periodic[0]:
        dx -= nx * np.round(dx / nx)
    dy = np.arange(ny) - Y
    if periodic[1]:
        dy -= ny * np.round(dy / ny)
    ox, oy = np.argsort(dx), np.argsort(dy)
    v = vtk["velocity"][zp][np.ix_(oy, ox)]
    fl = vtk["flag"][zp][np.ix_(oy, ox)]
    ux = np.where(fl == 0, v[..., 0], np.nan)
    uy = np.where(fl == 0, v[..., 1], np.nan)
    return dx[ox], dy[oy], dz, ux, uy


def oscillation_extrema(traj: dict, ny: int, R: float, U0: float, t_skip: float = 0.0, prominence=0.02):
    """Lateral extrema of a channel trajectory.

    Returns (T, yn, idx): T = t U0 / H, yn = (y - yc)/(H/2 - R) (so +-1 means
    touching a wall) and the indices of the local extrema of yn with T > t_skip.
    """
    from scipy.signal import find_peaks
    H = ny - 2
    yc = 0.5 * (ny - 1)
    T = traj["t"] * U0 / H
    yn = (traj["Y"] - yc) / (0.5 * H - R)
    pk, _ = find_peaks(yn, prominence=prominence)
    tr, _ = find_peaks(-yn, prominence=prominence)
    idx = np.sort(np.r_[pk, tr])
    return T, yn, idx[T[idx] > t_skip]



# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def plot_duct_validation(vtk, G, nu, convergence=None, title=None, path=None):
    """Profiles, 2D fields, error map and (optionally) a grid-convergence panel.

    convergence: optional (H_values, errors).
    """
    import matplotlib.pyplot as plt
    Y, Z, u, ue, err = duct_cross_section_error(vtk, G, nu)
    ny, nz = vtk["dims"][1], vtk["dims"][2]
    a, b = 0.5 * (ny - 2), 0.5 * (nz - 2)
    umax = np.nanmax(ue)
    zm, ym = nz // 2, ny // 2

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)
    s = np.linspace(-a, a, 400)
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(s / a, duct_velocity(s, Z[zm, 0], a, b, G, nu) / umax, "k-", label="analytical")
    ax.plot(Y[zm] / a, u[zm] / umax, "o", mfc="none", c="C3", ms=5, label="LBM")
    ax.set(xlabel="y / a", ylabel=r"$u_x/u_{max}$", title=f"Profile along y (z−z$_c$ = {Z[zm,0]:+.1f})")
    ax.legend(); ax.grid(alpha=.3)
    s = np.linspace(-b, b, 400)
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(s / b, duct_velocity(Y[0, ym], s, a, b, G, nu) / umax, "k-", label="analytical")
    ax.plot(Z[:, ym] / b, u[:, ym] / umax, "s", mfc="none", c="C0", ms=5, label="LBM")
    ax.set(xlabel="z / b", ylabel=r"$u_x/u_{max}$", title=f"Profile along z (y−y$_c$ = {Y[0,ym]:+.1f})")
    ax.legend(); ax.grid(alpha=.3)

    ax = fig.add_subplot(gs[0, 2])
    if convergence is not None:
        H, e = map(np.asarray, convergence)
        ax.loglog(H, e, "o-", c="C2", label="LBM (relative L2)")
        ax.loglog(H, e[0] * (H[0] / H) ** 2, "k--", lw=1, label="slope −2")
        ax.set_xticks(H); ax.set_xticklabels([str(int(h)) for h in H]); ax.minorticks_off()
        ax.set(xlabel="channel width H", ylabel="relative L2 error", title="Grid convergence")
        ax.legend(); ax.grid(alpha=.3, which="both")
    else:
        ax.axis("off")
        ax.text(0.1, 0.5, f"relative L2 error = {err:.2e}", fontsize=13)

    ext = [-a, a, -b, b]
    for k, (F, t) in enumerate([(u, "LBM"), (ue, "Analytical")]):
        ax = fig.add_subplot(gs[1, k])
        c = ax.imshow(F / umax, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=1)
        ax.contour(Y, Z, np.nan_to_num(F) / umax, levels=np.linspace(.1, .9, 9), colors="w", linewidths=.6)
        ax.set(title=rf"{t}  $u_x/u_{{max}}$", xlabel="y", ylabel="z"); fig.colorbar(c, ax=ax, shrink=.85)
    ax = fig.add_subplot(gs[1, 2])
    e = (u - ue) / umax
    m = np.nanmax(np.abs(e))
    c = ax.imshow(e, origin="lower", extent=ext, cmap="RdBu_r", vmin=-m, vmax=m)
    ax.set(title=r"(LBM − exact) / $u_{max}$", xlabel="y", ylabel="z"); fig.colorbar(c, ax=ax, shrink=.85)
    fig.suptitle(title or f"LBM vs analytical duct flow   (relative L2 error {err:.2e})", fontsize=13)
    if path:
        fig.savefig(path, dpi=120, bbox_inches="tight")
    return fig, err



def plot_squirmer_model(cases, R, B1, half_width=14.0, path=None):
    """Analytic Blake flow (top) vs LBM flow (middle) for several squirmers,
    slip profiles and on-axis velocity (bottom).

    cases: list of dicts with keys 'beta', 'label', 'vtk', 'state'.
    Returns (fig, errors) where errors[label] = {"axis": ..., "near_field": ...}:
    relative L2 differences LBM vs Blake of u_x along the swimming axis
    (R < |x| < 3R) and of (u_x, u_y) over the plane for 1.25R < r < 3R.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    U0 = 2 * B1 / 3
    L = half_width
    g = np.linspace(-L, L, 241)
    GX, GY = np.meshgrid(g, g)
    th = np.linspace(0, 2 * np.pi, 17)[:-1]
    colors = ["C3", "k", "C0", "C2", "C1"]
    n = len(cases)
    fig = plt.figure(figsize=(4.6 * n, 12.5))
    gs = fig.add_gridspec(3, n, height_ratios=[1, 1, .85], hspace=.32, wspace=.18)
    axis_data, axis_errors = [], {}
    for k, c in enumerate(cases):
        beta = c["beta"]
        ax = fig.add_subplot(gs[0, k])
        u = blake_velocity(GX, GY, np.zeros_like(GX), R, B1, beta)
        ax.pcolormesh(GX, GY, np.hypot(u[0], u[1]) / U0, cmap="magma_r", vmin=0, vmax=2.2, shading="auto")
        ax.streamplot(g, g, np.nan_to_num(u[0]), np.nan_to_num(u[1]), density=1.3, color="0.25",
                      linewidth=.6, arrowsize=.7)
        ax.add_patch(Circle((0, 0), R, fc="w", ec="k", lw=1.2, zorder=3))
        for t in th:
            rh = np.array([np.cos(t), np.sin(t), 0.0])
            us = B1 * (1 + beta * rh[0]) * (rh[0] * rh - np.array([1, 0, 0]))
            ax.arrow(R * rh[0], R * rh[1], 170 * us[0], 170 * us[1], head_width=.7, color="deepskyblue",
                     ec="k", lw=.4, length_includes_head=True, zorder=4)
        ax.arrow(-1.8, 0, 3.0, 0, head_width=.6, color="k", zorder=5)
        ax.set(aspect="equal", xlim=(-L, L), ylim=(-L, L), xticks=[], yticks=[],
               title=f"{c['label']}\nanalytical (Blake), lab frame")

        ax = fig.add_subplot(gs[1, k])
        px, py, dz, ux, uy = squirmer_plane(c["vtk"], c["state"])
        m = np.abs(px) <= L + 1
        ax.pcolormesh(px[m], py, np.hypot(ux, uy)[:, m] / U0, cmap="magma_r", vmin=0, vmax=2.2, shading="nearest")
        ax.streamplot(np.linspace(px[m][0], px[m][-1], m.sum()), py, np.nan_to_num(ux[:, m]),
                      np.nan_to_num(uy[:, m]), density=1.3, color="0.25", linewidth=.6, arrowsize=.7)
        ax.add_patch(Circle((0, 0), R, fc="none", ec="k", lw=1.2, ls="--", zorder=3))
        ax.set(aspect="equal", xlim=(-L, L), ylim=(-L, L), xticks=[], yticks=[],
               title=f"LBM, plane z − z$_c$ = {dz:+.2f}")

        j = int(np.argmin(np.abs(py)))
        ua = blake_velocity(px, np.full_like(px, py[j]), np.full_like(px, dz), R, B1, beta)[0]
        sel = (np.abs(px) > R) & (np.abs(px) < 3 * R) & ~np.isnan(ux[j])
        PX, PY = np.meshgrid(px, py)
        ub = blake_velocity(PX, PY, np.full_like(PX, dz), R, B1, beta)
        rr = np.sqrt(PX ** 2 + PY ** 2 + dz ** 2)
        mk = (rr > 1.25 * R) & (rr < 3 * R) & ~np.isnan(ux)
        axis_errors[c["label"]] = {
            "axis": float(np.sqrt(np.sum((ux[j][sel] - ua[sel]) ** 2) / np.sum(ua[sel] ** 2))),
            "near_field": float(np.sqrt(np.sum((ux[mk] - ub[0][mk]) ** 2 + (uy[mk] - ub[1][mk]) ** 2)
                                        / np.sum(ub[0][mk] ** 2 + ub[1][mk] ** 2))),
        }
        axis_data.append((c, px, ux[j], py[j], dz, colors[k % len(colors)]))

    ax = fig.add_subplot(gs[2, 0])
    t = np.linspace(0, np.pi, 200)
    for c, *_, col in axis_data:
        ax.plot(np.degrees(t), np.sin(t) * (1 + c["beta"] * np.cos(t)), c=col, label=c["label"])
    ax.axhline(0, c=".5", lw=.6)
    ax.set(xlabel="θ from swimming direction [deg]", ylabel=r"$u_\theta/B_1$", title="Surface slip",
           xticks=[0, 45, 90, 135, 180]); ax.legend(fontsize=8); ax.grid(alpha=.3)

    ax = fig.add_subplot(gs[2, 1:])
    s = np.linspace(-20, 20, 400)
    for c, px, uxl, yoff, dz, col in axis_data:
        ua = blake_velocity(s, np.full_like(s, yoff), np.full_like(s, dz), R, B1, c["beta"])[0]
        ax.plot(s / R, ua / U0, c=col, lw=1.3, label=f"{c['label']} analytical")
        ok = ~np.isnan(uxl) & (np.abs(px) <= 20)
        ax.plot(px[ok] / R, uxl[ok] / U0, "o", ms=3.5, mfc="none", c=col, label=f"{c['label']} LBM")
    ax.axvspan(-1, 1, color=".85"); ax.axhline(0, c=".5", lw=.6)
    ax.set(xlabel="distance along swimming axis / R", ylabel=r"$u_x/U_0$", xlim=(-5, 5),
           title="Flow along the swimming axis: LBM vs analytical")
    ax.legend(fontsize=7.5, ncol=3, loc="lower center", bbox_to_anchor=(.5, -.45)); ax.grid(alpha=.3)
    fig.suptitle(r"Squirmer model: $\mathbf{u}_s = B_1(1+\beta\cos\theta)(\cos\theta\,\hat{\mathbf{r}}-\mathbf{e})$,"
                 r"  $U_0=2B_1/3$.  Colour: $|\mathbf{u}|/U_0$", fontsize=12.5)
    if path:
        fig.savefig(path, dpi=110, bbox_inches="tight")
    return fig, axis_errors


# def plot_trajectories(runs, ny, R, U0, title=None, path=None):
#     """Trajectories in the channel, lateral position and orientation vs time.

#     runs: list of (traj_dict, label, color).
#     """
#     import matplotlib.pyplot as plt
#     H = ny - 2
#     yc, ylo, yhi = 0.5 * (ny - 1), 0.5, ny - 1.5
#     fig = plt.figure(figsize=(13, 3.2 + 2.6 * (len(runs) + 2)))
#     gs = fig.add_gridspec(len(runs) + 2, 1, hspace=.55)
#     for k, (tr, lab, col) in enumerate(runs):
#         ax = fig.add_subplot(gs[k])
#         ax.axhspan(ylo - 3, ylo, color=".7"); ax.axhspan(yhi, yhi + 3, color=".7")
#         ax.axhline(yc, ls=":", c="k", lw=.8)
#         x = tr["X"] - tr["X"][0]
#         ax.plot(x, tr["Y"], c=col, lw=1.4)
#         j = np.linspace(0, len(x) - 1, 30).astype(int)
#         ax.quiver(x[j], tr["Y"][j], tr["ex"][j], tr["ey"][j], angles="xy", scale=40, width=.0025, color="k")
#         ax.set(ylim=(ylo - 2, yhi + 2), xlim=(-10, x[-1] + 10), ylabel="y",
#                title=f"{lab}: trajectory (y stretched; arrows = swimming direction)")
#     ax1 = fig.add_subplot(gs[-2]); ax2 = fig.add_subplot(gs[-1], sharex=ax1)
#     for tr, lab, col in runs:
#         T = tr["t"] * U0 / H
#         ax1.plot(T, (tr["Y"] - yc) / (0.5 * H - R), c=col, label=lab)
#         ax2.plot(T, tr["angle"], c=col, label=lab)
#     for s in (-1, 1):
#         ax1.axhline(s, c=".5", ls="--", lw=.8)
#     ax1.set(ylabel=r"$(y-y_c)/(H/2-R)$", ylim=(-1.15, 1.15), title="Lateral position (±1 = touching a wall)")
#     ax2.axhline(0, c="k", lw=.6)
#     ax2.set(ylabel="orientation [deg]", xlabel=r"time $tU_0/H$", title="Swimming direction")
#     for ax in (ax1, ax2):
#         ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper right")
#     if title:
#         fig.suptitle(title, fontsize=12)
#     if path:
#         fig.savefig(path, dpi=120, bbox_inches="tight")
#     return fig




    #######
def plot_trajectories(runs, ny, R, U0, title=None, path=None):
    """Trajectories in the channel, lateral position and orientation vs time.

    runs: list of (traj_dict, label, color).
    """
    import matplotlib.pyplot as plt
    H = ny - 2
    yc, ylo, yhi = 0.5 * (ny - 1), 0.5, ny - 1.5
    fig = plt.figure(figsize=(13, 3.2 + 2.6 * (len(runs) + 2)))
    gs = fig.add_gridspec(len(runs) + 2, 1, hspace=.55)
    for k, (tr, lab, col) in enumerate(runs):
        ax = fig.add_subplot(gs[k])
        ax.axhspan(ylo - 3, ylo, color=".5"); ax.axhspan(yhi, yhi + 3, color=".7")
        ax.axhline(yc, ls=":", c="k", lw=.8)
        x = tr["X"] - tr["X"][0]
        ax.plot(x, tr["Y"], c=col, lw=3.4)
        j = np.linspace(0, len(x) - 1, 30).astype(int)
        ax.quiver(x[j], tr["Y"][j], tr["ex"][j], tr["ey"][j], angles="xy", scale=10, width=.0025, color="k")
        ax.set(ylim=(ylo - 2, yhi + 2), xlim=(-10, x[-1] + 10), ylabel="y",
               title=f"{lab}: trajectory (y stretched; arrows = swimming direction)")
    ax1 = fig.add_subplot(gs[-2]); ax2 = fig.add_subplot(gs[-1], sharex=ax1)
    for tr, lab, col in runs:
        T = tr["t"] * U0 / H
        ax1.plot(T, (tr["Y"] - yc) / (0.5 * H - R), c=col, label=lab)
        ax2.plot(T, tr["angle"], c=col, label=lab)
    for s in (-1, 1):
        ax1.axhline(s, c=".5", ls="--", lw=.8)
    ax1.set(ylabel=r"$(y-y_c)/(H/2-R)$", ylim=(-1.15, 1.15), title="Lateral position (±1 = touching a wall)")
    ax2.axhline(0, c="k", lw=.6)
    ax2.set(ylabel="orientation [deg]", xlabel=r"time $tU_0/H$", title="Swimming direction")
    for ax in (ax1, ax2):
        ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper right")
    if title:
        fig.suptitle(title, fontsize=12)
    if path:
        fig.savefig(path, dpi=120, bbox_inches="tight")
    return fig

