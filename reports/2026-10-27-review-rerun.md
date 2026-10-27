# 2026-10-27: Rerun of all studies with the 2026-10-26 setup (after the paper review)

Branch `claude/review-rerun`. Setup (`experiments/methods_review.jl`): HyperCloud 0.5.0 / **HyperCloudCore 0.7.3** /
PDEStudioCore 0.4.2; SSP time steppers at every order (RK2, RK3 SSP; **RK4SSP** = SSPRK(5,4) of Spiteri–Ruuth for
K ≥ 4), the **AlignedRusanov** flux (AR; the plain Rusanov flux RU as control), the MLS (divergence) order frozen
at 2 (the review's candidate default), frozen at K, or coupled to the MOOD order (0); MOOD U2 with δ = 1e-4 Δx^D
instead of the strict δ = 0 (review point 4). Everything on the cluster's CPU nodes in Float64 (the CUDA path has
a single-precision effect, reports/2026-10-22-…). Method labels `<stepper>MUSCL<K>(<variant>,<AR|RU>,MLS<m>)`.
Queues and job ids: `cluster/queues/2026-10-26-*.txt` and `cluster/queues/logs/`; data (gitignored, flat
save-path layout) in `data/review/{convergence,mass,euler,vortex}/`.

Config changes from the review (all in this branch):
- **Mach3Step ghosts** (point 5): the ~60 misclassified ghosts were points on a boundary line up to rounding that
  fell through the tolerance bands of the tag function into the "fallback" outflow tag; the tag function is now
  exactly complementary to the interior test (0 mismatches at 600×200, checked with a tag-vs-nearest-face script).
- **ClainShock exact boundary states** (point 5): the four plane shocks between the quadrant states, with the
  Rankine–Hugoniot speeds −0.663 (TL/TR, BR/TR) and −0.422 (BL/BR, BL/TL), prescribe the ghost states on all four
  sides (`ExactQuadrantBC`); `clainshock_outflow.jl` keeps the extrapolating boundaries as control.
- **RK2Upwind1 → RK2MUSCL1** (K = 1) in all defaults.
- **Isentropic vortex** (point 3): new `experiments/Euler/IsentropicVortex.jl` (periodic [0, 10]², β = 5, one period
  t = 10, jitter 0.2, R = 2.5, MOOD on), error against the initial condition (`tools/vortex_error.jl`), Float64 and
  Float32.
- **Sod** (point 1): in the conservation study with the exact boundary fluxes (0, (1 − 0.1) t, 0), at R = 5.5 and 2.5.

## 1. Convergence (point 3: the range of the shock runs)

Three experiments (1-D advection N ≤ 3163, 2-D advection and Burgers N ≤ 1000), Float64, CFL 0.2, no limiting,
SSP steppers + AR, lattice and jitter 0.2, **R = 5.5 (as before) and R = 2.5 (the shock runs)**, K ≤ 4 at 2.5
(Remark 3.8). Mean EOC of the last three refinements (L1 error at N_max); full tables in
`reports/2026-10-27-convergence-review-tables.md`.

| Experiment | K | lattice, R 5.5 | lattice, R 2.5 | jitter 0.2, R 5.5 | jitter 0.2, R 2.5 |
|---|---|---|---|---|---|
| 1-D linear | 2 | 2.00 (1.1e-03 at 3163) | 2.00 (2.1e-04 at 3163) | 2.00 (1.0e-03 at 3163) | 2.00 (2.0e-04 at 3163) |
| 1-D linear | 3 | 3.00 (1.6e-05 at 3163) | 2.00 (9.9e-05 at 3163) | 3.00 (1.6e-05 at 3163) | 2.00 (9.8e-05 at 3163) |
| 1-D linear | 4 | 4.00 (2.1e-08 at 3163) | 4.00 (4.9e-10 at 3163) | 3.98 (2.1e-08 at 3163) | 2.98 (2.5e-08 at 3163) |
| 1-D linear | 5 | 5.00 (5.3e-10 at 3163) | – | 4.98 (5.4e-10 at 3163) | – |
| 2-D burgers | 2 | 2.15 (8.2e-04 at 1000) | 2.06 (1.6e-04 at 1000) | 2.14 (8.3e-04 at 1000) | 1.97 (2.0e-04 at 1000) |
| 2-D burgers | 3 | 2.93 (8.7e-05 at 1000) | 2.98 (7.7e-06 at 1000) | 2.93 (8.8e-05 at 1000) | 2.91 (8.5e-06 at 1000) |
| 2-D burgers | 4 | 4.30 (1.3e-06 at 1000) | 4.12 (4.2e-08 at 1000) | 4.28 (1.4e-06 at 1000) | 3.89 (5.7e-05 at 200) |
| 2-D burgers | 5 | 4.85 (5.3e-07 at 1000) | – | 4.85 (5.4e-07 at 1000) | – |
| 2-D linear | 2 | 2.08 (2.7e-03 at 1000) | 2.02 (5.4e-04 at 1000) | 2.08 (2.6e-03 at 1000) | 1.99 (5.4e-04 at 1000) |
| 2-D linear | 3 | 3.00 (1.9e-04 at 1000) | 3.00 (1.7e-05 at 1000) | 3.00 (1.9e-04 at 1000) | 2.91 (1.9e-05 at 1000) |
| 2-D linear | 4 | 4.03 (6.0e-07 at 1000) | 4.01 (2.0e-08 at 1000) | 4.02 (5.9e-07 at 1000) | 3.24 (8.2e-08 at 1000) |
| 2-D linear | 5 | 5.00 (7.5e-08 at 1000) | – | 4.99 (7.5e-08 at 1000) | – |

- Everything at R = 5.5 reproduces the earlier study (the stepper and flux change nothing on smooth data): 2 / 3 /
  4 / 5 on lattice and jittered points.
- **R = 2.5 on the lattice keeps the orders** (2.0 / 3.0 / 4.0–4.1 in 2-D) with **errors 5–30× smaller** than at
  R = 5.5 (e.g. 2-D advection MUSCL4: 2.0e-8 vs 6.0e-7 at 1000²), as the smaller stencil weights the near
  neighbours. In 1-D, MUSCL3 is only second order at R = 2.5 (both lattice and jitter): the degree-3 divergence
  fit needs 3.5 Δx in 1-D (docs), so it falls back to degree 2.
- **R = 2.5 on jittered points costs order for K = 4**: 2.98 in 1-D, 3.24 on 2-D advection (3.78 on Burgers, 1.1e-7 at 1000²), while K = 2, 3 keep 2.0 / 2.9. With 0.2 jitter some
  stencils of radius 2.5 Δx have too few points for the degree-3 fit and drop a degree. So the shock-run range is
  fine for K ≤ 3 and for K = 4 on lattices; on perturbed clouds K = 4 wants R ≈ 3.5–5.5 (or more particles in the
  neighbourhood).

## 2. Mass conservation (point 1)

141 CPU jobs (2333 runs) of the same matrix as reports/2026-10-24-mass-conservation.md, now with the SSP steppers
and the aligned flux (AR) and with Sod added. Conservation defect = (m(t_end) − m(0) − exact boundary flux)/m,
sampling-free (`tools/mass_table.jl`); full tables `reports/2026-10-27-mass-review-tables.md`.

**1-D Burgers shock (jitter 0.4, R = 5.5), defect at N = 3163 — the MLS-order question:**

| K | coupled MLS 0 + EPD1 | coupled MLS 0 + EPD2 | frozen MLS = K + EPD1 | frozen MLS = K + EPD2 | **frozen MLS 2 + EPD1** | **frozen MLS 2 + EPD2** | plain (no MOOD) |
|---|---|---|---|---|---|---|---|
| 2 | −3.9e-4 (= frozen 2) | −6.7e-4 | −3.9e-4 | −6.7e-4 | −3.9e-4 | −6.7e-4 | +1.1e-4 |
| 3 | **+0.16** (was +2.8e-2) | **+5.8e-2** (was +3.6e-2) | −8.1e-5 (was +5.4e-3) | +8.7e-4 (was +1.7e-3) | −4.3e-4 | −6.3e-4 | +2.3e-4 |
| 4 | **+0.16** (was +2.7e-2) | **+6.1e-2** (was +2.8e-2) | −1.5e-3 (was +8.4e-5) | +6.1e-4 (was +7.6e-4) | −4.6e-4 | −5.9e-4 | +5.1e-4 |
| 5 | **+1.0 (fails)** (was −1.4e-3) | **+0.15** (was −4.1e-3) | −7.0e-3 (was −1.2e-2) | −3.6e-4 (was −3.5e-4) | −4.9e-4 | −6.7e-4 | −2.3e-3 |

("was": the plain Rusanov / classical RK4 runs of 2026-10-24; the floor of this test is ±3e-4, see that report.)

- **The coupled MLS order must not be combined with the aligned flux**: its flat defect grows from 3 % to 16 %
  (EPD1) and 6 % (EPD2), and MUSCL5 + EPD1 loses all its mass. That is consistent with what the flux is for
  (the HyperCloudCore docstring: it assumes a *uniform* divergence degree with MOOD acting on the reconstruction
  only); on mixed-degree pairs the weight-sign orientation makes the mismatch systematic. The coupled mode is
  out, as the review asks (point 1).
- **Frozen MLS 2 (the review's candidate default) is conservative to the floor at every order** (−4 … −7e-4,
  the same as MUSCL2 itself), with EPD1 and EPD2 alike.
- Frozen MLS = K with the aligned flux: K = 3 is now at the floor with EPD1 (−8e-5, from +5e-3), K = 5 improves
  (−7e-3 from −1.2e-2) but stays flat, K = 4 with EPD1 gets slightly worse (−1.5e-3 from +8e-5). So freezing at K
  still leaves 1e-3 … 7e-3 for K ≥ 4; the mechanism (fallback particles with a high-degree divergence next to
  full-order ones) is only partly cured by the flux.
- Limiters and plain high order unchanged (VK/minmod 4e-4 … 4e-3 flat, plain MUSCL5 −2.3e-3): no drops, no
  difference from the flux.
- Jitter, range, α (MUSCL4 EPD1, frozen 4): same picture as before, except that **R = 3.5 no longer blows up**
  (flat +3e-3 instead of −100 % at N ≥ 2000); R = 4.5 −1.2e-3, 5.5 −1.5e-3, 7.5 −8e-4, 9.5 −5e-4; jitter 0.1 … 0.5:
  +2e-4 … −1.9e-3 (plain MUSCL4 +2e-4 … +8e-4, growing with jitter and not converging, as before).
- Smooth periodic 1-D/2-D: identical to 2026-10-24 (1-D 1e-8 … 1e-15; 2-D ~1e-4 at 200² with jitter, independent
  of order and limiter, second order; lattice exact). No MOOD drops, so the flux makes no difference.

**Sod (1-D Euler, jitter 0.2, N = 100 … 1996), defects of (mass, momentum, energy) at N = 1996 with the exact
boundary fluxes (0, 0.9 t, 0):**

| Method | R = 5.5: mass / momentum / energy | R = 2.5: mass / momentum / energy |
|---|---|---|
| RK2MUSCL1 (K = 1) | −7.6e-7 / −3.7e-4 / +8.8e-7 | +6.8e-5 / −9.5e-4 / +7.8e-5 |
| MUSCL2 EPD1, MLS 2 | −3.4e-6 / −3.0e-5 / −5.9e-6 | +6.2e-5 / −8.1e-5 / +6.4e-5 |
| MUSCL2 EPD2, MLS 2 | −8.5e-6 / −3.2e-5 / −1.0e-5 | +4.9e-5 / −1.0e-4 / +4.5e-5 |
| MUSCL3 EPD1, MLS 2 | −4.5e-6 / −2.6e-5 / −9.0e-6 | +6.8e-5 / −4.4e-5 / +6.8e-5 |
| MUSCL3 EPD2, MLS 2 | −4.7e-6 / −7.9e-6 / −6.7e-6 | +5.3e-5 / −9.5e-5 / +4.2e-5 |
| MUSCL3 EPD1, MLS 3 | −2.8e-7 / +4.7e-5 / +5.5e-6 | +1.4e-4 / **+2.1e-2** / +1.5e-4 |
| MUSCL3 EPD2, MLS 0 | +5.7e-3 / +4.2e-2 / +1.3e-2 | **fails (−100 %)** |
| MUSCL4 EPD1, MLS 2 | −8.4e-6 / −3.5e-5 / −1.4e-5 | +5.7e-5 / −8.1e-5 / +5.0e-5 |
| MUSCL4 EPD2, MLS 2 | −1.6e-6 / −3.1e-6 / −2.8e-6 | +7.2e-5 / +5.2e-6 / +7.9e-5 |
| MUSCL4 EPD1, MLS 4 | +1.9e-5 / +2.1e-4 / +3.6e-5 | +3.9e-4 / +3.4e-3 / +1.0e-3 |
| MUSCL4 EPD2, MLS 0 | +1.3e-3 / +1.7e-2 / +4.6e-3 | **fails (−100 %)** |

- At R = 5.5 the frozen-2 variants conserve mass, momentum and energy to 1e-6 … 3e-5 **and converge** (orders 2–4.5
  over N); frozen K leaves 2e-5 … 2e-4 (momentum worst); the coupled mode is flat at 1e-3 … 4e-2.
- At **R = 2.5 the coupled mode fails outright** (MUSCL3/4 + EPD2 + MLS 0 lose the whole mass), frozen K = 3/4
  leave a 2 % momentum defect (frozen MLS 3 at K = 3) resp. 3e-3, and even the frozen-2 variants and MUSCL1 sit on
  a flat 5 … 8e-5 (mass) that does not converge — a boundary effect of the 2.5 Δx ghost layer (the same floor for
  K = 1 … 4), not of the schemes. So Sod needs R ≥ 3.5 for a clean conservation statement, which answers the
  review's "why 5.5" for Sod.

## 3. Euler shock problems (points 1, 2, 5)

Mach3Step 900×300, DMR 1440×450 (coupled K = 4 at 960×300), ClainShock 800² with the exact boundaries, CFL 0.8,
Float64, CPU. Δmass = final mass relative to RK2MUSCL2(EPD1,AR,MLS2). Full tables:
`reports/2026-10-27-euler-review-tables.md` (`tools/euler_table.jl`).

| Method | Mach3Step: min p / Δmass / steps | DMR: min p / Δmass | ClainShock: min p / Δmass |
|---|---|---|---|
| MUSCL1 (K = 1) | 1.00 / +0.2 % / 13166 | 1.00 / −0.2 % | 0.029 / −0.0 % |
| MUSCL2 EPD1 / EPD2 / VK (AR, MLS 2) | 0.19 / 0.27 / 0.23, ±0.3 % | 1.00, ±0.0 % | 0.026–0.029, ±0.0 % |
| **MUSCL3 EPD1 / EPD2, frozen MLS 2** | 0.19 / 0.28, −0.1 / +0.3 % / 13.6k | 1.00, ±0.0 % | 0.021 / 0.029, ±0.0 % |
| MUSCL3 EPD1 / EPD2, frozen MLS 3 | 0.21 / 0.17, **−1.7 / −2.0 %** / 40.5k | 0.98 / 1.00, −0.1 / 0.0 % | 0.018 / 0.025, ±0.0 % |
| MUSCL3 EPD1 / EPD2, coupled MLS 0 | **fails (t ≈ 0.003)** | **fails (t ≈ 0.005)** | **fails (t ≈ 0.001)** |
| **MUSCL4 EPD1 / EPD2, frozen MLS 2** | 0.068 / 0.26, −0.2 / +0.1 % / 13.6k | 1.00, ±0.0 % | 0.027 / 0.029, ±0.0 % |
| MUSCL4 EPD1 / EPD2, frozen MLS 4 | 0.14 / (EPD2 not run to the end), **−2.1 %** / 41.3k | 1.00 / 0.99, 0.0 % | 0.022 / 0.028, ±0.0 % |
| MUSCL4 EPD1 / EPD2, coupled MLS 0 | **fails** | **fails** (960×300) | **fails** |
| MUSCL2/3/4 EPD2, MLS 2, **plain Rusanov** | = AR to 3 digits | = AR | = AR |

- **Positivity with a frozen divergence order is fixed** (review point 2): every frozen-K run (MLS 3 for K = 3,
  MLS 4 for K = 4) now completes with positive pressure on all three problems, where all of them blew up within
  ~500 steps with the old flux (REPORT_2026-10-21, …-22). The cost: 3× the steps (the CFL uses the degree-3/4
  divergence) and, on Mach3Step, a −1.7 … −2.1 % mass defect (on DMR and ClainShock ≤ 0.1 %).
- **Frozen MLS 2 + EPD2 — the review's candidate default — runs on the step and the DMR** (point 1) with positive
  pressure, mass within ±0.3 % of the second-order reference, and MUSCL2 cost (13.6k steps on the step).
- **The coupled mode fails everywhere with the aligned flux** (within 14 … 300 steps; the DMR K = 4 EPD2 at
  t = 0.145): it should not be offered with this flux.
- Aligned and plain Rusanov give identical results with frozen MLS 2 (the flux only differs on pairs with an order
  drop below the full divergence degree, which frozen 2 never has for K ≥ 2).
- ClainShock exact vs extrapolating boundaries (EPD2, MLS 2): final mass +0.3850 vs +0.3788 (K = 2), +0.3850 vs
  +0.3787 (K = 3), +0.3850 vs +0.3732 (K = 4): **with the exact boundary states the mass is the same for every
  order (to 4 digits), with the outflow boundaries it depends on the order (up to 3 %)** — the boundary was the
  error source the paper suspected (§7.5), now a result.
- Not run to the end: Mach3Step RK4SSPMUSCL4(EPD2,AR,MLS4) (would have needed ~5.8 h; skipped, rerun at a smaller
  grid if needed). Runtimes on 32 EPYC cores: MUSCL2 10–20 min, frozen-2 MUSCL3/4 25–85 min, frozen-K 1.5–5 h.

## 4. Isentropic vortex with MOOD (point 3)

Periodic [0, 10]², one period (t = 10), jitter 0.2, R = 2.5, CFL 0.8, N = 50 … 400; L1(ρ) error against the
initial condition (integral, ×100 of the mean), `reports/2026-10-27-vortex-review-table.md`.

| Method | Float64: error at 400² / order (last 3) | Float32: error at 400² / order |
|---|---|---|
| MUSCL3, no MOOD (MLS 3) | 3.0e-3 / 2.9 | 8.6e-3 / floor from 283² |
| MUSCL4, no MOOD (MLS 4) | 1.6e-5 / 3.6 | **floor 2–5e-3 from 141²** (single precision) |
| MUSCL2 EPD2 (MLS 2) | 0.53 / 1.0 | 0.52 / 1.0 |
| MUSCL3 EPD2 (MLS 2) | 0.62 / 0.9 | 0.63 / 0.9 |
| MUSCL4 EPD2 (MLS 2) | 0.63 / 0.8 | 0.63 / 0.8 |
| MUSCL4 EPD2 (MLS 4) | 0.44 / 1.0 | 0.48 / 0.8 |

- **MOOD does not leave smooth flow at full order** (review point 3): with EPD2 (U2, δ = 1e-4 Δx²) every scheme is
  **first order** with errors 100–40 000× larger than without MOOD; the vortex core is smeared (min ρ 0.58 at 400²
  vs the exact 0.49 that the unlimited schemes keep). MOOD flags the smooth extremum of the vortex; the U2
  smooth-extremum test does not rescue it here. This is the most important new result for the paper: the
  high-order claim for shock runs with MOOD needs either a working smooth-extremum detection (order maps would
  show where it drops — point 4) or a larger δ.
- Without MOOD the schemes reach 2.9 (K = 3) and 3.4–3.9 (K = 4) at R = 2.5 in Float64; **in Float32 the K = 4
  error saturates at 2–5e-3 from 141²** (the single-precision effect the review asked about in point 4), and
  with MOOD on Float32 and Float64 are indistinguishable.
- Mass: 1e-5 … 4e-4 with MOOD (jittered points), 1e-8 … 1e-10 without.

## 5. Summary and verdicts

| Review point | Result | Verdict |
|---|---|---|
| 1 coupled mode / frozen k_div = 2 | coupled MLS 0 + aligned flux: 6–16 % flat mass defect (1-D), fails on all 2-D problems and on Sod at R 2.5; frozen 2 + EPD2 runs on step and DMR, conservative to ±0.3 % | coupled mode not usable; **frozen 2 + EPD2 recommended** |
| 1 measure conservation directly | Sod m/momentum/E defects with exact fluxes; 2-D mass vs reference; vortex | done (§2, §3) |
| 2 positivity with frozen degree | frozen MLS 3/4 now positive on all three 2-D problems (aligned flux) | PASS, at 3× cost and −2 % mass on the step |
| 3 κ = 2.5 high order | lattice: full order, 5–30× smaller errors; jitter 0.2: K = 4 drops to 3.2–3.8 | PASS K ≤ 3; K = 4 needs R > 2.5 on perturbed clouds |
| 3 MOOD on smooth flow | vortex with EPD2: first order for every K | **FAIL** — MOOD degrades smooth flow |
| 3 SSP RK4 | SSPRK(5,4) used for all K = 4 runs | done |
| 4 single precision | vortex K = 4 without MOOD: Float32 floor 2–5e-3 | confirmed |
| 5 Mach 3 step ghosts | 0 misclassified (was ~60) | fixed |
| 5 RK2Upwind1 | replaced by RK2MUSCL1 | done |
| 5 2-D Riemann boundaries | exact plane-shock states; mass then independent of K (outflow: up to 3 % K-dependence) | done, attribution confirmed |

Overall: **INCONCLUSIVE** for the paper's high-order claims with MOOD (the vortex result), PASS for positivity and
the frozen-2 default.
