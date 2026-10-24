# 2026-10-24: Mass conservation on irregular point clouds (1-D Burgers shock, 1-D/2-D smooth periodic)

103 CPU jobs (Float64, 32 or 16 cores), 1-D sweeps N = 32 … 3163, 2-D sweeps N = 32 … 200; HyperCloud 0.5.0 /
HyperCloudCore 0.7.0 / PDEStudioCore 0.4.1. Data: `data/mass/` (flat, README with the job ids), every run in
`reports/2026-10-24-mass-conservation-tables.md`, CSV `data/mass/mass_runs.csv` (`tools/mass_table.jl`).
Job script `cluster/run_mass.sh` (jitter, weight range R, weight α, MLS order, method, N range per job).

## 0. How mass is measured

The hyperbolic stat preset stores the discrete mass m(t) = Σ_i V_i u_i (`mass`) and `relative_mass` = m / mass
of the analytical solution *sampled on the particles*. For a discontinuous reference the sampling has an O(Δx)
error (2e-4 at N = 3163 on the shock problem, the same for every scheme), which would hide the conservation
defect of the scheme. The tables therefore use the **conservation defect** computed from m(t) alone,

    defect = (m(t_end) − m(0) − exact boundary flux) / m(t_end),

with the exact flux 0 for periodic domains and (f(u_L) − f(u_R)) t_end = 0.5 t_end for the outflow Riemann
problem (shock stays inside). A conservative scheme gives 0 up to rounding. (`relative_mass` is kept in the full
tables for comparison; for the smooth periodic cases both agree.)

Problems: **1-D Burgers shock** (`conservation_burgers`: u_L = 1 → u_R = 0 at x = −2, outflow BCs, t = 10, jitter
0.4, R = 5.5, MLS order fixed = scheme order, MOOD U2 where used), **1-D linear advection** of a Gaussian, periodic
(`Convergence1D`, t = 10), **2-D Burgers** of a Gaussian before shock formation, periodic (`ConvergenceBurgers`,
t = 0.5). "Order" below is the order of |defect| over the last three refinements (negative = growing).

## 1. The dominant effect: MOOD order drops with a coupled MLS order (MLS 0)

1-D shock, jitter 0.4, R = 5.5, defect at N = 3163 (and its behaviour over N):

| Method | MLS order | defect | over N |
|---|---|---|---|
| RK4MUSCL3(EPD1) | **0** | **−2.8e-2** | flat 2.5–3.3e-2 from N = 80 on |
| RK4MUSCL3(EPD2) | **0** | **+3.6e-2** | flat |
| RK4MUSCL4(EPD1) | **0** | **+2.7e-2** | flat (2.7e-2 for every N ≥ 252) |
| RK4MUSCL4(EPD2) | **0** | **+2.8e-2** | flat |
| RK4MUSCL5(EPD1) | 0 | −1.4e-3 | flat |
| RK4MUSCL5(EPD2) | 0 | −4.1e-3 | flat |
| RK4MUSCL2(EPD1/EPD2) | 0 | −3.8e-4 / −5.5e-4 | = fixed MLS 2 (identical: at order 2 nothing can drop below the degree-1 divergence) |

A 3 % mass error that does not change with N at all: the fallback particles (order 1 after a MOOD drop) use a
different divergence stencil than their neighbours, the pairwise antisymmetry of the fluxes is lost at every
interface of the shock zone, and the number of such interfaces per shock is fixed. This is the mechanism behind
the +3 … +10 % of MLS 0 + EPD1 on the Euler problems (REPORT_2026-10-21.md, reports/2026-10-22-…). The sign
depends on the problem (gain here for MUSCL4, loss for MUSCL3).

## 2. Fixed MLS order with MOOD, limiters, no limiting

Same setting (jitter 0.4, R = 5.5), fixed MLS = scheme order, defect at N = 3163:

| | order 2 | order 3 | order 4 | order 5 |
|---|---|---|---|---|
| plain (no limiter, no MOOD) | +1.1e-4 | +2.3e-4 | +5.1e-4 | −2.3e-3 |
| EPD1 (MOOD U2) | −3.8e-4 | +5.4e-3 | +8.4e-5 | −1.2e-2 |
| EPD2 (MOOD U2) | −5.5e-4 | +1.7e-3 | +7.6e-4 | −3.5e-4 |
| VK | −4.3e-4 | +5.1e-4 | −7.2e-4 | −3.6e-3 |
| minmod | -3.60e-04 | +2.0e-3 | +1.5e-3 | −3.8e-3 |

- **Floor of the test:** on the regular lattice (jitter 0) the plain and EPD1 schemes of order 2 and 4 all give
  the *same* defect, −2.0e-4 at N = 3163, O(Δx) (order 1.00), and exactly 0 (1e-13) at some N (100, 159, 200,
  631–1259). So this is not interior non-conservation but a boundary/shock-position effect of the outflow setup
  that depends on how the shock sits on the lattice; the 1-D shock test cannot resolve defects below ~3e-4.
  Everything in the 1e-4 range in the tables is at that floor.
- Above the floor and **not converging**: all MOOD/limiter variants of order 3 and 5 (5e-4 … 1.2e-2), minmod at
  every order (1.5–3.8e-3), and the plain fifth-order scheme (−2.3e-3). With a fixed MLS order the drop to
  order 1 (or the limiter) changes the reconstruction but the divergence weights stay antisymmetric only if the
  *same* divergence degree is used on both sides, which is the case; the remaining defect must come from the
  neighbourhoods where the degree-5 (or 3) fit does not have enough points at jitter 0.4 and falls back (the
  docs: a fit with fewer neighbours than unknowns lowers its degree), again mixing stencils. The even orders
  2 and 4 with EPD1/VK stay at the floor.
- MUSCL5 with fixed MLS 5 and EPD1 loses 1.2 % mass (−9e-3 … −1.4e-2 for all N ≥ 252), the worst fixed-order
  case; EPD2 at order 5 is at 3.5e-4. Odd orders are the problematic ones here, as they were for positivity.

## 3. Irregularity (jitter)

1-D shock, R = 5.5, fixed MLS, defect at N = 3163:

| jitter | MUSCL2 plain | MUSCL2 EPD1 | MUSCL4 plain | MUSCL4 EPD1 |
|---|---|---|---|---|
| 0 | −2.0e-4 (floor) | −2.0e-4 | −2.0e-4 | −2.0e-4 |
| 0.1 | +1.9e-4 | +1.7e-4 | +2.2e-4 | +1.6e-4 |
| 0.2 | +1.5e-4 | +9.9e-5 | +2.6e-4 | +6.9e-7 |
| 0.3 | +1.0e-4 | −5.5e-5 | +3.1e-4 | +3.4e-5 |
| 0.4 | +1.1e-4 | −3.8e-4 | +5.1e-4 | +8.4e-5 |
| 0.5 | +2.2e-4 | −5.5e-4 | +8.0e-4 | −6.2e-4 |

The jitter alone costs at most a few 1e-4 on the shock problem (plain MUSCL4 grows from 2e-4 to 8e-4 with the
jitter and does *not* converge with N, orders −1 … −4), i.e. the same magnitude as the floor; the MOOD
variants sit on top of that.

**1-D smooth periodic advection** (exact mass constant), defect at N = 3163:

| | MUSCL2 | MUSCL3 | MUSCL4 | MUSCL5 |
|---|---|---|---|---|
| jitter 0.2 | +1.6e-8 (order 3.5) | −1.4e-10 (6.6) | +3.3e-13 (4.3) | −6e-16 |
| jitter 0.4 | +7.2e-8 (order 2.0) | −9.6e-10 (4.2) | +2.3e-12 (2.3) | +4.7e-14 (6.9) |

So for smooth solutions in 1-D the jitter-induced non-conservation converges with (at least) the scheme order
and is negligible (1e-8 … 1e-15); the Oct-17 "10 % mass loss on jittered grids" was the box (discontinuity +
MOOD), not the jitter.

**2-D smooth periodic Burgers** (exact mass constant), defect at N = 200²:

| jitter | R | MUSCL2 | MUSCL2(EPD1) | MUSCL2(VK) | MUSCL4 | MUSCL4(EPD1) | MUSCL4(VK) |
|---|---|---|---|---|---|---|---|
| 0 | 3.5 / 5.5 | 0 | 1e-16 | 0 | 0 | 1e-16 | 1e-16 |
| 0.2 | 3.5 | +8.3e-5 | +7.5e-5 | +7.9e-5 | +1.1e-4 | +1.1e-4 | +1.1e-4 |
| 0.2 | 5.5 | +1.4e-4 | +1.3e-4 | +1.3e-4 | +1.1e-4 | +1.1e-4 | +1.1e-4 |
| 0.4 | 3.5 | +2.2e-4 | +2.1e-4 | +2.2e-4 | +2.2e-4 | +2.2e-4 | +2.2e-4 |
| 0.4 | 5.5 | +2.6e-4 | +2.5e-4 | +2.6e-4 | +2.3e-4 | +2.3e-4 | +2.4e-4 |

Orders of the defect 1.4–2.3 for every variant (e.g. MUSCL4, jitter 0.4, R 5.5: −7.5e-3 at 32², −7.8e-4 at
100², +2.2e-4 at 200²). In 2-D the jitter alone *does* break conservation at the 1e-4 level at 200² with a
defect that is **independent of the scheme order and of the limiter** (MUSCL2 ≈ MUSCL4, plain ≈ EPD1 ≈ VK),
roughly proportional to the jitter (0.2 → 0.4 doubles it), only weakly dependent on R, and converging at about
second order: the asymmetry c_ij V_i ≠ −c_ji V_j of the MLS stencils on irregular 2-D points (REPORT_2026-10-17
§3). The lattice is exact to rounding. This is far larger than in 1-D (1e-12 for MUSCL4 at N = 3163).

## 4. Interpolation range and weight

1-D shock, jitter 0.4, EPD1, fixed MLS, defect at N = 3163:

| R | MUSCL2(EPD1) | MUSCL4(EPD1) |
|---|---|---|
| 3.5 | −6.6e-4 | **−1.0 (blows up for N ≥ 1996; −2e-2 before)** — too few neighbours for the degree-3 fit in 1-D |
| 4.5 | −5.0e-4 | −4.8e-3 (flat over N) |
| 5.5 | −3.8e-4 | +8.4e-5 |
| 7.5 | −4.7e-5 | +1.1e-3 (flat) |
| 9.5 | −1.7e-4 | +1.1e-3 (flat) |

| α (MUSCL4 EPD1, R 5.5) | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|
| defect | +1.8e-4 | +8.4e-5 | −6.3e-5 | −4.7e-4 |

For MUSCL4 the range has an optimum around 5.5: smaller ranges starve the degree-3/4 fits (R = 3.5 is unusable,
R = 4.5 leaves a flat 5e-3), larger ranges bring the defect back to a flat 1e-3 (more neighbours of different
MOOD order inside each stencil). The weight sharpness α matters little (all within the floor); a sharper
weight (α = 4) is slightly worse.

## 5. Summary

1. **Keep the MLS (divergence) order fixed.** MLS 0 with MOOD costs a flat 3 % of the mass on a single 1-D shock,
   independent of N; with a fixed order the even-order schemes (2, 4) with EPD1/EPD2/VK are conservative to the
   floor of the test (3e-4), the odd orders (3, 5) and minmod leave flat defects of 5e-4 … 1.2e-2.
2. **Jitter alone is harmless in 1-D** (1e-8 … 1e-15 on smooth data, ≤ 8e-4 and comparable to the test floor on
   the shock) **but costs ~1e-4 at 200² in 2-D** on smooth data, independent of order and limiter, converging at
   ~second order and proportional to the jitter. The lattice is exactly conservative.
3. **Range:** 5.5 Δx is the sweet spot for MUSCL4; 3.5 is unusable, 4.5 and ≥ 7.5 leave flat 1e-3 … 5e-3 defects
   with MOOD. α is a second-order effect.
4. The 1-D shock test itself has an O(Δx) floor from the outflow/shock-position setup (−2e-4 at N = 3163 on the
   lattice, zero at some N); a conservation study of the *interior* scheme is cleaner on periodic problems.

## 6. Open points

- RK4MUSCL2(minmod), 1-D shock: the minmod limiter at order 2 is very slow (30 min per run at N ≥ 2000); the
  job hit its 2 h limit, the missing N = 1996, 2512, 3163 were rerun as single jobs (included).
- Why the lattice shock runs give exactly 0 at some N and O(Δx) at others (boundary flux vs shock position) was
  not traced; it does not affect the comparisons above (all at the same N and setup).
- Not covered: Euler (the Mach3Step/DMR/ClainShock mass results are in the earlier reports), moving grids,
  other weight functions.
