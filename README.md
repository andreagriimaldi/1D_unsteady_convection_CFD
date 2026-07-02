# Inviscid Burgers Equation — 1D Finite-Volume Solver

Explicit finite-volume solution of the inviscid Burgers equation in conservative form,

```
∂u/∂t + ∂/∂x (u²/2) = 0,   −L ≤ x ≤ L,   t ≥ 0
```

using first-order upwind (donor-cell) fluxes and an adaptive time step at Courant number `C = 1`. Final exam, *Introducción a la Dinámica de Fluidos Computacional* (FIE-UNDEF, SiCaNLab).

**Author:** Andrea Grimaldi,
**Course instructor:** Prof. Edgardo A. Serafín

---

## Method

The conservative form is discretized on a cell-centred mesh of `N` control volumes with two ghost cells for the Dirichlet boundaries. Each interior cell is advanced explicitly:

```
uᵢⁿ⁺¹ = uᵢⁿ − (Δt/Δx) (Fₑ − F_w),   F = u²/2
```

Face fluxes use upwind (donor-cell) selection keyed to the sign of the local wave speed `f'(u) = u`, not the sign of the flux `F` (which is always ≥ 0 and direction-blind). For the non-negative fields of both configurations the donor is always the western cell.

The time step is recomputed every iteration from the current field,

```
Δt = C · Δx / max|u|,   C = 1,
```

and clipped so the march lands exactly on the requested output times. `C = 1` sits at the CFL stability limit and, for this scheme, is also the point of minimum numerical diffusion.

## Configurations

| | Configuration 1 | Configuration 2 |
|---|---|---|
| Initial condition | Gaussian `u₀ = exp(−2x²)` | Square wave: `u = 2` on `[−2.9, −2]`, else `1` |
| Domain | `L = 2` | `L = 3` |
| Boundaries | `u_L = u_R = 0` | `u_L = u_R = 1` |
| Behaviour | shock born from smooth data at `t_b ≈ 0.824` | shock (speed `3/2`) + trailing rarefaction fan |

## Results

**Configuration 1 — Gaussian profile.** The forward face steepens into a shock while the rear face opens into a rarefaction; the peak is convected downstream and decays after the shock forms.

![Configuration 1: u(x) at selected times](case1.png)

**Configuration 2 — Square wave.** The leading edge propagates as a shock at the Rankine–Hugoniot speed `s = 3/2`; the trailing edge spreads into a rarefaction fan that erodes the plateau after the interaction time `t = 1.8`.

![Configuration 2: u(x) at selected times](case2.png)

Full analysis — characteristic and Riemann verification, numerical-diffusion study, mesh independence, and first-order shock-position convergence — is in [`Grimaldi_final_report.pdf`](Grimaldi_final_report.pdf).

## Repository contents

```
.
├── 1D_convection.py          # solver (both configurations, plots, optional GIF)
├── Report.pdf                # technical report
├── case1.png                 # Config 1: u(x) at selected times
├── case1_umax.png            # Config 1: max(u(x)) over time
├── case1_xA.png              # Config 1: u(xA) over time
├── case2.png                 # Config 2: u(x) at selected times
├── case2_umax.png            # Config 2: max(u(x)) over time
├── case2_xA.png              # Config 2: u(xA) over time
├── .gitignore
└── README.md
```

## Usage

Requirements: Python 3 with NumPy and Matplotlib.

```bash
pip install numpy matplotlib
python 1D_convection.py
```

Select the case at the top of the script:

```python
config = 1          # 1 → Gaussian, 2 → square wave
make_animation = 0  # 1 → export burgers_case{config}.gif
nvol = 5000         # number of control volumes
```

Running produces three figures per configuration: `u(x)` at the output times, `max u` versus time, and the probe history `u(x_A)` versus time (`x_A = 1.2`).

## Notes

- **Conservative form is mandatory.** Only the flux form reproduces the correct shock speed; the non-conservative form `u ∂u/∂x` is undefined across a discontinuity.
- **Outflow boundary is inert.** For rightward flow the last interior cell draws its flux from the west, so the eastern ghost value `u_R` is never read; any front exits cleanly without reflection.
- **First-order limitation.** Numerical diffusion `Γ_num ~ |u|Δx/2` smears fronts on coarse meshes; it vanishes linearly in `Δx` and is controlled by mesh refinement (`N = 5000` at convergence).