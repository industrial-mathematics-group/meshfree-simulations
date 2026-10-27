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
- **R = 2.5 on jittered points costs order for K = 4**: 2.98 in 1-D, 3.24 on 2-D advection (3.9 on Burgers up to
  200²; the N ≥ 251 line of that sweep is still running), while K = 2, 3 keep 2.0 / 2.9. With 0.2 jitter some
  stencils of radius 2.5 Δx have too few points for the degree-3 fit and drop a degree. So the shock-run range is
  fine for K ≤ 3 and for K = 4 on lattices; on perturbed clouds K = 4 wants R ≈ 3.5–5.5 (or more particles in the
  neighbourhood).

## 2. Mass conservation (point 1)

MASS_SECTION

## 3. Euler shock problems (points 1, 2, 5)

EULER_SECTION

## 4. Isentropic vortex with MOOD (point 3)

VORTEX_SECTION

## 5. Summary and verdicts

SUMMARY_SECTION
