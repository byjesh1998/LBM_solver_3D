
# lbm-squirmer 

A compact 3D **lattice Boltzmann (LBM)** solver for simulating a spherical microswimmer (**squirmer**) swimming in a rectangular channel.

The code uses a D3Q19 lattice, TRT collision, no-slip walls, and a fully resolved spherical swimmer. It can be used to study **pushers, pullers, and neutral swimmers** near walls.

<p align="center">
  <img src="docs/figures/pusher.png" width="85%"><br>
  <img src="docs/figures/puller.png" width="85%"><br>
  <img src="docs/figures/neutral.png" width="85%"><br>
  <em>Trajectories of different squirmer between the walls of a slit channel; Orientations of squirmer are colour coded.</em>
</p>


## `Features`

- D3Q19 lattice, two-relaxation-time (TRT) collision with the "magic" parameter Λ = 3/16
- Guo body forcing (pressure-driven background flow)
- No-slip walls in y and/or z (halfway bounce-back); x is always periodic
- Squirmer with the Blake/Lighthill two-mode slip velocity (parameter β selects pusher, neutral or puller)
- Weak, smooth wall potential with a hard safety gap
- Validated: duct flow to second-order accuracy, bulk swimming speed and near-field flow vs. Blake's solution
- OpenMP parallel; binary VTK output for ParaView


---

## Contents
1. [Start running the simulation](#1-start-running-the-simulation)
2. [Repository layout](#2-repository-layout)
3. [Theory and numerical method](#3-theory-and-numerical-method)
4. [Input file reference](#4-input-file-reference)
5. [Output files](#5-output-files)
6. [Validation](#6-validation)
7. [Tests](#7-tests)
8. [Tutorial notebook](#8-tutorial-notebook)
9. [Units and dimensionless numbers](#9-units-and-dimensionless-numbers)
10. [Limitations and known issues](#10-limitations-and-known-issues)
11. [References](#11-references)

---

## 1. Start running the simulation

### Requirements

* `C++17 compiler` `GCC ≥ 9 or Clang ≥ 10` `OpenMP (optional)`

For tests and post-processing:

* `Python ≥ 3.9` `NumPy` `SciPy` `Matplotlib` `pytest` `Jupyter`

### Build

```bash
> git clone https://github.com/byjesh1998/LBM_solver_3D.git
> cd LBM_solver_3D

> make
```

Or with CMake:

```bash
>cmake -B build
>cmake --build build
```

This creates:

```text
bin/lbm_squirmer
```

---

## Run a simulation

### Duct-flow validation

```bash
> ./bin/lbm_squirmer inputs/duct.in
```

### Neutral squirmer

```bash
> ./bin/lbm_squirmer inputs/channel_neutral.in
```

### Puller with custom parameters

```bash
> ./bin/lbm_squirmer inputs/channel_puller.in beta=2 angle=15 prefix=my_run
```

Parameters can be changed directly from the command line: `key=value`

For example:

```bash
> ./bin/lbm_squirmer inputs/channel_neutral.in beta=-3 radius=6
```

### Run with OpenMP

```bash
> export OMP_NUM_THREADS=8
> ./bin/lbm_squirmer inputs/channel_neutral.in
```

---


## 2. Repository layout

```
LBM_solver_3D/
├── src/
│   ├── main.cpp              driver: read input, time loop, output
│   └── utils/
│       ├── lattice.hpp       D3Q19 velocities, weights, equilibrium, node types
│       ├── params.hpp/.cpp   parameters, input-file parser
│       ├── lbm.hpp/.cpp      TRT collision + streaming + bounce-back, squirmer coupling
│       ├── squirmer.hpp/.cpp squirmer model and parameters, rigid-body dynamics, wall potential
│       ├── io.hpp/.cpp       VTK, trajectory and state output
│       └── analytic.hpp/.cpp exact duct / plane Poiseuille solutions
├── inputs/                   ready-to-run input files (*.in)
├── outputs/                  simulation results (standard test results are included)
├── postprocessing/
│   └── lbm_post.py           run the solver, read output, exact solutions, plots
├── tests/                    pytest suite (duct flow, bulk squirmer, trajectories)
├── notebooks/tutorial.ipynb  step-by-step tutorial
├── examples/reference_data/  long-run trajectories used by the notebook
├── docs/figures/             figures used in this README
├── Makefile, CMakeLists.txt
└── README.md
```

---

## 3. Theory and numerical method  

All quantities described in simulations are in **lattice units**: $\Delta x = \Delta t = 1$, reference density $\rho_0 = 1$.

### $\color{purple}{\text{3.1 Lattice Boltzmann equation}}$

The fluid is described by particle distribution functions $f_i(\mathbf x, t)$ moving with discrete
velocities $\mathbf c_i$, $i = 0,\dots,18$ (D3Q19):

$$
f_i(\mathbf x + \mathbf c_i, t+1) = f_i(\mathbf x, t) + \Omega_i(\mathbf x,t) + S_i(\mathbf x,t).
$$

where $\Omega_i$ is the collision operator and $S_i$ external force  term. The D3Q19 consider 19 discrete velocity points in 3D space, and each velocities $c_{i}$ weighed with weights $w_{i}$, as described in following table: 
<table border="0">
<tr>
<td>
  
| group | $\mathbf c_i$ | weight $w_i$ |
|---|---|---|
| rest | $(0,0,0)$ | $1/3$ |
| 6 face neighbours | $(\pm1,0,0), (0,\pm1,0), (0,0,\pm1)$ | $1/18$ |
| 12 edge neighbours | $(\pm1,\pm1,0), (\pm1,0,\pm1), (0,\pm1,\pm1)$ | $1/36$ |

 </td>
<td>

 <img src="docs/figures/d3q19.png" width="75%"><br>

</td>
</tr>
</table>


The speed of sound is $c_s^2 = 1/3$. Density and momentum are moments of distribution function $f_i$
(with a half-force correction, see 3.3):

$$
\rho = \sum_i f_i,\qquad \rho\mathbf u = \sum_i f_i\,\mathbf c_i + \tfrac12\mathbf F .
$$

The equilibrium distribution is chosen as the second-order expansion of the Maxwell–Boltzmann distribution:

$$
f_i^{eq} = w_i\,\rho\left[1 + \frac{\mathbf c_i\cdot\mathbf u}{c_s^2} + \frac{(\mathbf c_i\cdot\mathbf u)^2}{2c_s^4} - \frac{\mathbf u\cdot\mathbf u}{2c_s^2}\right],
$$

which come from the local equilibrium approximation and is necessary to recover Navier-Stokes equation from Boltzmann equation. In the low-Mach-number limit this recovers the incompressible Navier–Stokes equations:

$$
\nabla\cdot\mathbf u = 0, \qquad
\rho\left(\partial_t\mathbf u + \mathbf u\cdot\nabla\mathbf u\right) = -\nabla p + \mu\nabla^2\mathbf u + \mathbf F,
$$

with pressure $p = c_s^2\rho$ and kinematic viscosity

$$
\nu = c_s^2\left(\tau - \tfrac12\right) = \frac{\tau - 1/2}{3}.
$$

### $\color{purple}{\text{3.2 Two-relaxation-time (TRT) collision}}$ 

Each population is split into symmetric and antisymmetric parts with respect to its opposite
direction $\bar{i}$ ($\mathbf{c}_{\bar{i}} = -\mathbf{c}_i$):

$$
f_i^\pm = \tfrac12\left(f_i \pm f_{\bar{i}}\right),\qquad
f_i^{eq,+} = w_i\rho\left[1 + \tfrac92(\mathbf c_i\cdot\mathbf u)^2 - \tfrac32 u^2\right],\qquad
f_i^{eq,-} = 3w_i\rho\,\mathbf c_i\cdot\mathbf u ,
$$

and each part relaxes at its own rate:

$$
\Omega_i = -\omega^+ \left(f_i ^+ - f_i^{eq,+}\right) - \omega^- \left(f_i ^- - f_i ^{eq,-}\right),\qquad
\omega^+ = \frac1\tau .
$$

$\omega^-$ is fixed by the "magic" parameter

$$
\Lambda = \left(\frac{1}{\omega^+} - \frac12\right)\left(\frac{1}{\omega^-} - \frac12\right) = \frac{3}{16}.
$$

With $\Lambda = 3/16$ the no-slip plane of halfway bounce-back sits exactly half-way between nodes,
independent of $\tau$. This makes TRT exact for plane Poiseuille flow and removes the
viscosity-dependent wall slip of the single-relaxation-time (BGK) model.

### $\color{purple}{\text{3.3  Body force (Guo forcing)}}$ 

A uniform force density $\mathbf F$ (e.g. a pressure gradient $-\mathrm dp/\mathrm dx$) enters as a source term,
split into symmetric and antisymmetric parts:

$$
S_i = \left(1-\tfrac{\omega^+}{2}\right)S_i^+ + \left(1-\tfrac{\omega^-}{2}\right)S_i^-,\qquad
S_i^+ = w_i\left[9(\mathbf c_i\cdot\mathbf u)(\mathbf c_i\cdot\mathbf F) - 3\,\mathbf u\cdot\mathbf F\right],\qquad
S_i^- = 3w_i\,\mathbf c_i\cdot\mathbf F .
$$

### $\color{purple}{\text{3.4 Walls: halfway bounce-back}}$ 

The outermost node layer in y (and optionally z) is solid. A population that would stream from a
fluid node $\mathbf x_f$ into a solid node is reflected back:

$$
f_{\bar i}(\mathbf x_f, t+1) = f_i^*(\mathbf x_f, t) - 2 w_i\,\rho\,\frac{\mathbf c_i\cdot\mathbf u_b}{c_s^2}.
$$

Here $f_i^*$ is the post-collision value and $\mathbf u_b$ the velocity of the boundary at the link midpoint
$\mathbf x_b = \mathbf x_f + \tfrac12\mathbf c_i$ (zero for the walls).
The no-slip plane lies at $\mathbf x_b$, so with fluid nodes $y = 1,\dots,N_y-2$
the walls are at $y = 0.5$ and $y = N_y - 1.5$ and the channel width is

$$
H = N_y - 2 .
$$

### $\color{purple}{\text{3.5 Squirmer model}}$ 

The squirmer model is a simplified model of swimming microswimmers in fluid. A squirmer is a rigid sphere of radius $R$ with orientation $\mathbf e$ whose surface drives the
fluid with a prescribed tangential slip velocity (Lighthill 1952; Blake 1971). Keeping the first two
modes, with $\cos\theta = \mathbf e\cdot\hat{\mathbf r}$, gives the surface slip velocity as follows:

$$
\mathbf u_s(\theta) = \left(B_1\sin\theta + B_2\sin\theta\cos\theta\right)\hat{\boldsymbol\theta}
= B_1\left(1 + \beta\cos\theta\right)\left(\cos\theta\,\hat{\mathbf r} - \mathbf e\right),
\qquad \beta = \frac{B_2}{B_1}.
$$

| β | swimmer | far field | propulsion |
|---|---|---|---|
| β < 0 | **pusher** (e.g. *E. coli*) | fluid expelled along the axis, drawn in at the sides | thrust from behind |
| β = 0 | **neutral** (e.g. *Volvox*) | source dipole, decays as $r^{-3}$ | — |
| β > 0 | **puller** (e.g. *Chlamydomonas*) | fluid drawn in along the axis, expelled at the sides | thrust from the front |

In unbounded Stokes flow the swimming speed is $\mathbf{U}_0 = \tfrac23 B_1 \mathbf e$, independent of β,
and the lab-frame flow is (Blake 1971)

$$
\mathbf u(\mathbf r) =
B_1\frac{R^3}{r^3}\left[(\mathbf e\cdot\hat{\mathbf r})\hat{\mathbf r} - \frac{\mathbf e}{3}\right]
+B_2\left(\frac{R^4}{r^4} - \frac{R^2}{r^2}\right)\frac{3(\mathbf e\cdot\hat{\mathbf r})^2 - 1}{2}\,\hat{\mathbf r}
+B_2\frac{R^4}{r^4}(\mathbf e\cdot\hat{\mathbf r})\left[(\mathbf e\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf e\right].
$$

The $B_2$ term is a force dipole (stresslet) decaying as $r^{-2}$, which dominates hydrodynamic
interactions with walls.

<p align="center"><img src="docs/figures/squirmer_model.png" width="90%"><br>
<em>Top: analytic flow fields around squirmer. Middle: LBM flow fields. Bottom: slip profile and on-axis velocity.</em></p>

### $\color{purple}{\text{3.6 Fluid–particle coupling}}$  

**Boundary condition.** Lattice nodes inside the sphere are flagged `PARTICLE`. Links from fluid
nodes into the sphere use the bounce-back rule of §3.4 with the local surface velocity

$$
\mathbf u_b = \mathbf U + \boldsymbol\Omega\times\mathbf r_b + \mathbf u_s(\hat{\mathbf r}_b),
\qquad \mathbf r_b = \mathbf x_f + \tfrac12\mathbf c_i - \mathbf X ,
$$

where $\mathbf X$ is the centre (minimum-image convention in periodic directions).

**Hydrodynamic force and torque.** These are obtained by Galilean-invariant momentum exchange
(Wen et al. 2014), summed over all boundary links:

$$
\mathbf F_h = \sum_{\text{links}}\left[\left(f_i^* + f_{\bar i}\right)\mathbf c_i - \left(f_i^* - f_{\bar i}\right)\mathbf u_b\right],
\qquad
\mathbf T_h = \sum_{\text{links}} \mathbf r_b\times\left[\cdots\right].
$$

**Covering and uncovering (Aidun et al. 1998).** When the sphere moves onto a fluid node, the
node's momentum $\sum_i f_i\mathbf c_i$ is given to the particle. When it uncovers a node, the node is
filled with $f_i^{eq}(\bar\rho, \mathbf U + \boldsymbol\Omega\times\mathbf r)$ ($\bar\rho$: mean density of
the fluid neighbours), and that momentum is removed from the particle.

**Rigid-body dynamics.** Explicit Newton–Euler integration each time step:

$$
\mathbf U \leftarrow \mathbf U + \frac{\mathbf F_h + \mathbf F_c + \mathbf F_w}{M},\quad
\boldsymbol\Omega \leftarrow \boldsymbol\Omega + \frac{\mathbf T_h + \mathbf T_c}{I},\quad
\mathbf X \leftarrow \mathbf X + \mathbf U,\quad
\mathbf e \leftarrow \frac{\mathbf e + \boldsymbol\Omega\times\mathbf e}{\lVert\mathbf e + \boldsymbol\Omega\times\mathbf e\rVert},
$$

with $M = \tfrac43\pi R^3\rho_p$ and $I = \tfrac25 M R^2$. With `planar = 1` the motion is restricted
to the x-y plane ($U_z = 0$, $\boldsymbol\Omega = \Omega_z\hat{\mathbf z}$). By symmetry this is an
exact solution of the equations of motion, but it can be unstable in 3D (see §10).

**Mass correction.** Moving-boundary bounce-back on a staircase surface does not conserve mass
exactly. Every `mass_correction_every` steps the deficit is redistributed isotropically,
$f_i \leftarrow f_i + w_i(1 - \langle\rho\rangle)$, over all fluid nodes.

### 3.7 Soft wall potential

The gap between the sphere surface and the wall is only one or two lattice spacings when the
swimmer is close, which the fluid solver cannot resolve (lubrication). A weak, smooth potential
acting along the wall normal prevents unphysical contact. For a surface-to-wall gap $h$:

$$
V(h) = \frac{\varepsilon F_S h_c}{3}\left(1 - \frac{h}{h_c}\right)^3,\qquad
F_w(h) = -V'(h) = \varepsilon F_S\left(1 - \frac{h}{h_c}\right)^2 \quad (h < h_c),
$$

and zero for $h \ge h_c$. The strength is measured in units of the Stokes drag at the swimming
speed, $F_S = 6\pi\mu R U_0$, so ε = 0.2 means at most 20 % of that drag. The force is central, so it
exerts no torque. As a last resort, a hard gap `hmin` is enforced: if $h < h_{min}$ the sphere is
placed back at $h = h_{min}$ and its wall-normal velocity is removed. The number of such
steps is reported as `contacts`.

### 3.8 Algorithm (one time step)

1. For every fluid node: compute ρ, **u**; TRT collision with forcing; push-stream each population;
   on links into solids apply bounce-back and accumulate momentum exchange.
2. Soft wall force → Newton–Euler update of **U**, **Ω**, **X**, **e** → hard gap.
3. Re-map the sphere on the lattice (covering/uncovering with momentum corrections).
4. Every `mass_correction_every` steps: global mass correction.

---

## 4. Input file reference

Input files contain `key = value` lines; `#` starts a comment. The same keys are accepted on the
command line (`key=value`, overriding the file). Booleans accept `0/1`, `true/false`, `yes/no`, `on/off`.

| key | default | description |
|---|---|---|
| `nx`, `ny`, `nz` | 64, 34, 32 | grid size (includes the wall layers) |
| `tau` | 1.0 | relaxation time, ν = (τ − ½)/3; must be > 0.5 |
| `lambda` | 0.1875 | TRT magic parameter Λ |
| `fx`, `fy`, `fz` | 0 | body-force density (background pressure-driven flow) |
| `walls_y` | 1 | no-slip walls at y = 0.5 and y = ny − 1.5 |
| `walls_z` | 0 | no-slip walls at z = 0.5 and z = nz − 1.5 |
| `steps` | 10000 | number of time steps |
| `tol` | 0 | stop when the relative change of Σu_x per `print_every` < tol (no squirmer only) |
| `print_every` | 1000 | console output interval |
| `traj_every` | 50 | trajectory output interval |
| `vtk_every` | 0 | VTK snapshot interval (0: final field only) |
| `mass_correction_every` | 10 | mass-correction interval (0: off) |
| `output_dir`, `prefix` | `outputs`, `run` | output files are `output_dir/prefix_*` |
| `squirmer` | 0 | enable the squirmer |
| `radius` | 4 | squirmer radius R (≥ 4 recommended) |
| `b1` | 0.015 | first squirming mode B₁; U₀ = 2B₁/3 (keep U₀ ≲ 0.02) |
| `beta` | 0 | β = B₂/B₁ (< 0 pusher, > 0 puller) |
| `rho_p` | 1 | particle/fluid density ratio |
| `x0`, `y0`, `z0` | domain centre | initial position (−1 means centre) |
| `angle` | 30 | initial angle from the x-axis in the x-y plane, **must satisfy \|angle\| < 45°** |
| `planar` | 1 | restrict motion to the x-y plane |
| `eps` | 0.2 | wall-potential strength (max force / F_S) |
| `hc` | 2.0 | wall-potential range (lattice units) |
| `hmin` | 0.5 | hard minimum surface-wall gap |


---

## 5. Output files

All files are written to `output_dir` with the chosen `prefix`.

| file | content |
|---|---|
| `prefix_final.vtk`, `prefix_<step>.vtk` | legacy binary VTK (structured points): `flag` (0 fluid, 1 wall, 2 particle), `density`, `velocity` |
| `prefix_traj.dat` | columns `t X Y Z ex ey ez Ux Uy Uz Wx Wy Wz Fwall_y Fwall_z`; X is unwrapped along x |
| `prefix_squirmer.dat` | squirmer state at the end of the run (`key value` pairs) |
| `prefix_summary.txt` | diagnostics: mean density, wall force, total body force, L2 error, contacts |
| `prefix_profile_y.dat`, `prefix_profile_z.dat` | centre-line profiles `coordinate u_lbm u_exact` (no squirmer) |

Python readers are in `postprocessing/lbm_post.py`: `read_vtk`, `read_traj`, `read_keyvalue`.

---

## 6. Validation

### 6.1 Duct flow (no squirmer)

Flow in a square duct driven by a body force $G$ is compared with the exact series solution for
$-a<y<a$, $-b<z<b$ (e.g. White, *Viscous Fluid Flow*):

$$
u_x(y,z) = \frac{16 a^2 G}{\mu\pi^3}\sum_{n=1,3,5,\dots}(-1)^{\frac{n-1}{2}}
\left[1 - \frac{\cosh(n\pi z/2a)}{\cosh(n\pi b/2a)}\right]\frac{\cos(n\pi y/2a)}{n^3},
\qquad a = b = \frac{H}{2}.
$$

<p align="center"><img src="docs/figures/duct_validation.png" width="85%"></p>

| channel width H | relative L2 error |
|---|---|
| 8 | 4.0 × 10⁻³ |
| 16 | 8.8 × 10⁻⁴ |
| 32 | 2.1 × 10⁻⁴ |

The error converges at second order. At steady state the force on the walls equals the total body
force to machine precision, and plane Poiseuille flow (walls in y only) is reproduced to round-off (L2 error 8 × 10⁻¹³, TRT with Λ = 3/16).

### 6.2 Squirmer in a periodic box

With R = 4 and U₀ = 0.01 (Re = 0.48) in a 40³–48³ periodic box:

- A neutral squirmer swims at 0.985 U₀, in a straight line with no rotation (machine precision).
- At this Re, pushers (β = −3) swim at 1.14 U₀ and pullers (β = +3) at 0.86 U₀, consistent with the known effect of fluid inertia.
- The near-field flow agrees with the Blake solution to 15–25 % (relative L2 over 1.25R < r < 3R).
  The remaining difference comes from the coarse sphere, periodic images and finite Re.

### 6.3 Squirmer in a slit channel (H = 32, R = 4, 2R/H = 0.25, Re = 0.48)

| swimmer | behaviour | notes |
|---|---|---|
| neutral, β = 0 | sustained wall-to-wall oscillation | period ≈ 7.4 H/U₀, wavelength ≈ 6.9 H, tilt ±14° |
| puller, β = +1 | damped oscillation towards the centreline | decay rate ≈ 0.02 U₀/H, period ≈ 14 H/U₀ |
| pusher, β = −3 | wall-to-wall motion, pressed against the walls | relies on the hard gap; sensitive to the wall model |

<p align="center"><img src="docs/figures/channel_trajectories.png" width="85%"><br>
<em>Neutral, puller and pusher trajectories (reference data, see <code>examples/reference_data/</code>).</em></p>

<p align="center"><img src="docs/figures/puller_beta1_channel.png" width="85%"><br>
<em>Puller β = +1: exponential fit of the decaying oscillation.</em></p>

---

## 7. Tests

```bash
python -m pytest tests -v             # fast: duct flow + bulk squirmer (~1.5 min)
python -m pytest tests -v --runslow   # also the channel-trajectory tests (~15 min)
make test / make test-all             # same, via make
```

| file | what is tested |
|---|---|
| `tests/test_duct_flow.py` | duct L2 error < 5×10⁻⁴; second-order convergence; wall force = body force; mass; exact plane Poiseuille; figure |
| `tests/test_squirmer_bulk.py` | neutral swimming speed; straight swimming; pusher > neutral > puller speed; mass; near-field flow vs Blake; figure |
| `tests/test_trajectories.py` (slow) | neutral squirmer visits both walls; puller oscillation decays; figure |

Figures are written to `tests/output/`.

---

## 8. Tutorial notebook

`notebooks/tutorial.ipynb` walks through the method step by step:
building and running the code, duct flow validation, the squirmer model, bulk validation and channel
trajectories. Short runs are executed live; long-time trajectories are loaded from
`examples/reference_data/` (set `RUN_LONG = True` in the notebook to recompute them).

---

## 9. Units and dimensionless numbers

The physics is set by dimensionless groups; lattice values are chosen to keep the Mach number
small and the sphere resolved.

| quantity | definition | typical value here |
|---|---|---|
| swimming Reynolds number | $Re = 2RU_0/\nu$ | 0.48 |
| confinement | $\kappa = 2R/H$ | 0.25 |
| squirmer parameter | $\beta = B_2/B_1$ | −3 … +3 |
| Mach number | $Ma = U_0/c_s = \sqrt3\,U_0$ | 0.017 |
| time unit | $H/U_0$ | 3200 steps |

To convert to physical units, choose a length scale $\Delta x$ (e.g. $R_{phys}/R$) and match the
viscosity, $\Delta t = \nu_{lat}\,\Delta x^2/\nu_{phys}$. Real microswimmers have $Re \sim 10^{-4}$–$10^{-2}$;
lower $Re$ is reached by increasing τ or decreasing $B_1$, at the cost of longer runs.

---

## 10. Limitations and known issues

- **Near-wall resolution.** When a swimmer is within ~1 lattice unit of a wall, lubrication is not
  resolved; outcomes near walls (especially for pushers) depend on `eps`, `hc`, `hmin` and on R.
  Check key results with a larger radius (R = 6–8) or with an explicit lubrication correction.
- **Out-of-plane instability.** Without `planar = 1`, a swimmer started in the x-y plane may drift
  in z because small lattice asymmetries grow.
- **Finite inertia.** Re ≈ 0.5 by default. Behaviour at Re → 0 may differ quantitatively.
- **Explicit coupling.** Newton–Euler integration is explicit; very light particles
  (ρ_p/ρ_f ≪ 1) can become unstable.
- **Single squirmer**, no restart files, single-node (OpenMP) parallelism.

---

## 11. References

- M. J. Lighthill, *On the squirming motion of nearly spherical deformable bodies through liquids at very small Reynolds numbers*, Commun. Pure Appl. Math. **5**, 109 (1952).
- J. R. Blake, *A spherical envelope approach to ciliary propulsion*, J. Fluid Mech. **46**, 199 (1971).
- T. Ishikawa, M. P. Simmonds, T. J. Pedley, *Hydrodynamic interaction of two swimming model micro-organisms*, J. Fluid Mech. **568**, 119 (2006).
- A. J. C. Ladd, *Numerical simulations of particulate suspensions via a discretized Boltzmann equation*, J. Fluid Mech. **271**, 285 (1994).
- C. K. Aidun, Y. Lu, E.-J. Ding, *Direct analysis of particulate suspensions with inertia using the discrete Boltzmann equation*, J. Fluid Mech. **373**, 287 (1998).
- Z. Guo, C. Zheng, B. Shi, *Discrete lattice effects on the forcing term in the lattice Boltzmann method*, Phys. Rev. E **65**, 046308 (2002).
- I. Ginzburg, F. Verhaeghe, D. d'Humières, *Two-relaxation-time lattice Boltzmann scheme*, Commun. Comput. Phys. **3**, 427 (2008).
- B. Wen, C. Zhang, Y. Tu, C. Wang, H. Fang, *Galilean invariant fluid–solid interfacial dynamics in lattice Boltzmann simulations*, J. Comput. Phys. **266**, 161 (2014).
- L. Zhu, E. Lauga, L. Brandt, *Low-Reynolds-number swimming in a capillary tube*, J. Fluid Mech. **726**, 285 (2013).
- T. Krüger et al., *The Lattice Boltzmann Method: Principles and Practice*, Springer (2017).
- F. M. White, *Viscous Fluid Flow*, McGraw-Hill.

---





             ┌───────────────-┐
             │ f distributions│
             └───────┬───────-┘
                     ↓
             calculate rho, u
                     ↓
                 collision
                     ↓
                 streaming
                     ↓
           ┌─────────┴─────────┐
           ↓                   ↓
        fluid cell          solid cell
           ↓                   ↓
      move normally       bounce back
                               ↓
                         force/torque






                         LBM
                          │
       ┌──────────────────┼──────────────────┐
       │                  │                  │
     FLUID             BOUNDARIES         PARTICLE
       │                  │                  │
       ↓                  ↓                  ↓
   f, rho, u          flag, wallForce       sq
       │                  │                  │
       └──────────────┬───┴──────────────────┘
                      ↓
                    step()
                      ↓
              fluid-particle coupling





# lbm-squirmer

A compact 3D **lattice Boltzmann (LBM)** solver for simulating a spherical **squirmer** swimming in a channel.

The code uses a D3Q19 lattice, TRT collision, no-slip walls, and a fully resolved spherical swimmer. It can be used to study **pushers, pullers, and neutral swimmers** near walls.

<p align="center">
  <img src="docs/figures/pusher.png" width="80%">
  <br>
  <img src="docs/figures/puller.png" width="80%">
  <br>
  <img src="docs/figures/neutral.png" width="80%">
  <br>
  <em>Example trajectories of pusher, puller, and neutral squirmers.</em>
</p>

---

## Features

* 3D **D3Q19 lattice Boltzmann** method
* **Two-relaxation-time (TRT)** collision
* Guo forcing for pressure-driven flow
* No-slip walls using halfway bounce-back
* Fully resolved spherical squirmer
* Pusher, neutral, and puller swimmers
* Optional soft wall repulsion
* Moving particle with force and torque coupling
* OpenMP parallelisation
* VTK output for ParaView
* Python post-processing and plotting
* Automated validation and tests
* Jupyter tutorial

---

## Quick start

### Requirements

You need:

* C++17 compiler
* GCC ≥ 9 or Clang ≥ 10
* OpenMP (optional)

For tests and post-processing:

* Python ≥ 3.9
* NumPy
* SciPy
* Matplotlib
* pytest
* Jupyter

### Build

```bash
git clone https://github.com/byjesh1998/LBM_solver_3D.git
cd LBM_solver_3D

make
```

Or with CMake:

```bash
cmake -B build
cmake --build build
```

This creates:

```text
bin/lbm_squirmer
```

---

## Run a simulation

### Duct-flow validation

```bash
./bin/lbm_squirmer inputs/duct.in
```

### Neutral squirmer

```bash
./bin/lbm_squirmer inputs/channel_neutral.in
```

### Puller with custom parameters

```bash
./bin/lbm_squirmer inputs/channel_puller.in beta=2 angle=15 prefix=my_run
```

Parameters can be changed directly from the command line:

```bash
key=value
```

For example:

```bash
./bin/lbm_squirmer inputs/channel_neutral.in beta=-3 radius=6
```

### Run with OpenMP

```bash
export OMP_NUM_THREADS=8
./bin/lbm_squirmer inputs/channel_neutral.in
```

---

## What is a squirmer?

A squirmer is a simple model of a swimming microorganism.

The swimmer is a sphere with a prescribed tangential surface velocity. The parameter `beta` controls the type of swimmer:

| `beta` | Type        | Example         |
| -----: | ----------- | --------------- |
|  `< 0` | **Pusher**  | *E. coli*       |
|    `0` | **Neutral** | *Volvox*        |
|  `> 0` | **Puller**  | *Chlamydomonas* |

The swimming speed in the ideal unbounded Stokes-flow model is

$$
U_0 = \frac{2}{3}B_1.
$$

The main squirmer parameters are:

* `b1` — swimming mode \(B_1\)
* `beta` — \(B_2/B_1\)
* `radius` — swimmer radius
* `angle` — initial swimming direction

---

## Numerical method

The solver uses:

* D3Q19 lattice
* TRT collision
* Guo forcing
* Halfway bounce-back boundaries
* Momentum-exchange fluid-particle coupling
* Newton-Euler particle dynamics
* A smooth wall potential to prevent unresolved contact

All simulations use lattice units:

$$
\Delta x = \Delta t = 1,
\qquad
\rho_0 = 1.
$$

The kinematic viscosity is

$$
\nu = \frac{\tau - 1/2}{3}.
$$

The TRT magic parameter is set to

$$
\Lambda = \frac{3}{16}.
$$

This gives good wall placement for halfway bounce-back and is used for the duct-flow validation.

For the full derivation and implementation details, see the documentation in `docs/`.

---

## Typical parameters

Some important input parameters are:

| Parameter | Default | Description                  |
| --------- | ------: | ---------------------------- |
| `nx`      |      64 | Grid size in x               |
| `ny`      |      34 | Grid size in y               |
| `nz`      |      32 | Grid size in z               |
| `tau`     |     1.0 | LBM relaxation time          |
| `fx`      |       0 | Body force in x              |
| `walls_y` |       1 | Enable walls in y            |
| `walls_z` |       0 | Enable walls in z            |
| `steps`   |   10000 | Number of time steps         |
| `radius`  |       4 | Squirmer radius              |
| `b1`      |   0.015 | First squirmer mode          |
| `beta`    |       0 | Pusher/puller parameter      |
| `rho_p`   |       1 | Particle/fluid density ratio |
| `angle`   |      30 | Initial swimming angle       |
| `planar`  |       1 | Restrict motion to x-y       |
| `eps`     |     0.2 | Wall-potential strength      |
| `hc`      |     2.0 | Wall-potential range         |
| `hmin`    |     0.5 | Minimum wall gap             |

The complete parameter list is available in the input files and source code.

---

## Output

Simulation output is written to the selected output directory.

For example:

```text
outputs/
├── run_final.vtk
├── run_traj.dat
├── run_squirmer.dat
├── run_summary.txt
├── run_profile_y.dat
└── run_profile_z.dat
```

### VTK files

VTK files contain:

* `flag`
* `density`
* `velocity`

They can be opened directly with **ParaView**.

### Trajectory file

`*_traj.dat` contains:

```text
t X Y Z ex ey ez Ux Uy Uz Wx Wy Wz Fwall_y Fwall_z
```

The x-position is unwrapped, so long trajectories can be analysed without jumps caused by periodic boundaries.

---

## Repository structure

```text
LBM_solver_3D/
├── src/
│   ├── main.cpp
│   └── utils/
│       ├── lattice.hpp/.cpp
│       ├── params.hpp/.cpp
│       ├── lbm.hpp/.cpp
│       ├── squirmer.hpp/.cpp
│       ├── io.hpp/.cpp
│       └── analytic.hpp/.cpp
│
├── inputs/
│   └── *.in
│
├── postprocessing/
│   └── lbm_post.py
│
├── tests/
│   ├── test_duct_flow.py
│   ├── test_squirmer_bulk.py
│   └── test_trajectories.py
│
├── notebooks/
│   └── tutorial.ipynb
│
├── examples/
│   └── reference_data/
│
├── docs/
│   └── figures/
│
├── Makefile
├── CMakeLists.txt
└── README.md
```

---

## Validation

The solver includes several validation cases.

### Duct flow

The LBM velocity profile is compared with the analytical solution for flow through a square duct.

The measured relative L2 error is:

| Channel width |    Relative L2 error |
| ------------: | -------------------: |
|             8 | \(4.0\times10^{-3}\) |
|            16 | \(8.8\times10^{-4}\) |
|            32 | \(2.1\times10^{-4}\) |

The error shows approximately **second-order convergence**.

---

### Squirmer in a periodic box

For a squirmer with \(R=4\) and \(U_0=0.01\):

* Neutral squirmer: approximately \(0.985U_0\)
* Pusher, \(\beta=-3\): approximately \(1.14U_0\)
* Puller, \(\beta=+3\): approximately \(0.86U_0\)

The near-field flow agrees with the analytical Blake solution to approximately **15–25% relative L2 error** for the tested range.

The remaining difference is attributed to the coarse spherical representation, periodic images, and finite Reynolds number.

---

### Squirmer in a channel

For \(H=32\), \(R=4\), and \(Re=0.48\):

| Swimmer              | Observed behaviour               |
| -------------------- | -------------------------------- |
| Neutral, \(\beta=0\) | Wall-to-wall oscillation         |
| Puller, \(\beta=+1\) | Damped oscillation toward centre |
| Pusher, \(\beta=-3\) | Strong wall interaction          |

<p align="center">
  <img src="docs/figures/channel_trajectories.png" width="80%">
  <br>
  <em>Example channel trajectories.</em>
</p>

---

## Tests

Run the fast tests with:

```bash
python -m pytest tests -v
```

Run all tests, including the longer trajectory tests:

```bash
python -m pytest tests -v --runslow
```

Or:

```bash
make test
make test-all
```

The tests cover:

* Duct-flow accuracy
* Second-order convergence
* Wall force balance
* Mass conservation
* Plane Poiseuille flow
* Bulk squirmer swimming
* Pusher/puller swimming speeds
* Near-field flow
* Channel trajectories

Test figures are written to:

```text
tests/output/
```

---

## Tutorial

A step-by-step tutorial is provided in:

```text
notebooks/tutorial.ipynb
```

The notebook covers:

1. Building the solver
2. Running a simulation
3. Duct-flow validation
4. The squirmer model
5. Bulk validation
6. Channel trajectories
7. Plotting the results

Long trajectory data is provided in:

```text
examples/reference_data/
```

Set:

```python
RUN_LONG = True
```

to recompute the long simulations.

---

## Units and dimensionless numbers

The simulations use lattice units.

Some typical dimensionless parameters are:

| Quantity           | Definition        | Typical value |
| ------------------ | ----------------- | ------------: |
| Reynolds number    | \(Re=2RU_0/\nu\)  |          0.48 |
| Confinement        | \(\kappa=2R/H\)   |          0.25 |
| Squirmer parameter | \(\beta=B_2/B_1\) |      −3 to +3 |
| Mach number        | \(Ma=U_0/c_s\)    |         0.017 |

The default simulations therefore operate at finite, but relatively low, Reynolds number.

Real microswimmers generally operate at lower Reynolds numbers. Lower \(Re\) can be reached by reducing \(B_1\) or increasing \(\tau\), at the cost of longer simulations.

---

## Limitations

This is a research-oriented numerical model, so some limitations are important.

### Near-wall resolution

Lubrication effects are not fully resolved when the swimmer is very close to a wall.

Results near contact can depend on:

* `radius`
* `eps`
* `hc`
* `hmin`

For important results, repeat simulations with a larger swimmer radius or use an explicit lubrication correction.

### Finite Reynolds number

The default Reynolds number is approximately:

$$
Re \approx 0.5.
$$

Therefore, results may differ quantitatively from the \(Re\rightarrow0\) Stokes-flow limit.

### Out-of-plane motion

Without:

```text
planar = 1
```

small lattice asymmetries can cause a swimmer initially placed in the x-y plane to drift in z.

### Particle dynamics

The particle dynamics use explicit Newton-Euler integration. Very light particles may become unstable.

### Other limitations

Currently the code supports:

* One squirmer
* No restart files
* OpenMP parallelism on a single node

---

## References

The implementation is based on the following works:

* Lighthill, M. J. (1952). *On the squirming motion of nearly spherical deformable bodies through liquids at very small Reynolds numbers.*
* Blake, J. R. (1971). *A spherical envelope approach to ciliary propulsion.*
* Guo, Z., Zheng, C., & Shi, B. (2002). *Discrete lattice effects on the forcing term in the lattice Boltzmann method.*
* Ginzburg, I., Verhaeghe, F., & d'Humières, D. (2008). *Two-relaxation-time lattice Boltzmann scheme.*
* Aidun, C. K., Lu, Y., & Ding, E.-J. (1998). *Direct analysis of particulate suspensions with inertia using the discrete Boltzmann equation.*
* Wen, B. et al. (2014). *Galilean invariant fluid-solid interfacial dynamics in lattice Boltzmann simulations.*
* Krüger, T. et al. (2017). *The Lattice Boltzmann Method: Principles and Practice.*

See the source code and original README for the complete reference list.

---

## License

Add your project license here.

For example:

```text
MIT License
```

if the project is released under MIT.

---

## Summary

`lbm-squirmer` provides a compact way to simulate a spherical microswimmer in a 3D channel using the lattice Boltzmann method.

It includes:

**LBM → squirmer → walls → particle dynamics → validation → analysis**

The code is intended primarily for research, teaching, and experimentation with microswimmer–wall interactions.

              
