# 2026-10-30: 2-D Euler benchmarks (experiments/Benchmarks2D) on the GPU in single precision

HyperCloud 0.5.1 / HyperCloudCore 0.7.8. Tests chosen from `/workspace/2D-Euler-Benchmarks.pdf` (the minimal suite there:
vortex, 2-D Riemann problems, DMR, forward step, Quirk, implosion, Mach 2000 jet; DMR, the Mach 3 step and the vortex
are already in experiments/). New here: Lax–Liu configurations 4, 6, 12, 15, 17, the implosion, Noh, Sedov, the Mach 80
and Mach 2000 jets, the Gresho vortex (classic and Mach 0.01), Kelvin–Helmholtz and the shock–vortex interaction.
Not done: Quirk and the hypersonic cylinder (need special point sets / curved walls), Rayleigh–Taylor (needs gravity),
multi-component problems, the steady airfoil cases.

**Setup.** 14 configs × the 12 methods known to work on 0.7.8: MUSCL1; MUSCL2 with VK, EPD1, EPD2; MUSCL3/MUSCL4 with
EPD1/EPD2 and the MLS divergence degree frozen at 2 or at K (AlignedRusanov flux, MOOD U2, δ = 1e-4 Δx²). Jitter 0.2
(the implosion on a lattice), R = 2.5, CFL 0.4. **168 runs on the cluster GPUs in Float32** (queues
`cluster/queues/2026-10-29-bench2d-gpu*.txt`, `2026-10-30-bench2d-gpu-3.txt`, results lines in
`cluster/queues/logs/*bench2d*-results.txt`), every run that ended with negative pressure **rerun in Float64 on the CPU**
(18 runs, locally with `tools/local_queue.sh`, queue `cluster/queues/2026-10-30-bench2d-f64-local.txt`).
All data in `data/`; `@sim <config>` (and `<config>_Float64` for the reruns) finds every run (`tools/check_config.jl`).
Tables: `tools/bench_table.jl <configs>`.

Checks per run (`tools/run_sim.jl` summary): reaches t_end, no NaN, ρ > 0 and p > 0 everywhere at t_end, no blow-up.
Benchmark statistics (post-processing): Noh and Gresho errors against the exact solution, implosion x↔y symmetry error,
Sedov shock radius vs the self-similar prediction, Kelvin–Helmholtz amplitude of the seeded mode.

## Summary

| Test | What it stresses | Float32 GPU | Float64 rerun | Verdict |
|---|---|---|---|---|
| Lax–Liu 4, 6, 12, 15, 17 (800²) | shock/contact/rarefaction interaction | 60/60 ok | – | PASS |
| Implosion (400², t = 2.5) | x↔y symmetry, long run | 12/12 ok; symmetry error 4e-7 (MUSCL1), 3e-6 (VK), 1e-4 … 9e-4 (MOOD) | – | PASS (MOOD breaks the symmetry 100× more than VK) |
| Noh (400², γ = 5/3) | strong shock, cold gas, wall heating | 7/12 ok; 5 MOOD runs with p ≈ −1e-7 (round-off of the cold gas), VK p = −1.5e-2 | MOOD all ok; VK still p < 0 | PASS for MOOD (needs Float64), FAIL for VK |
| Sedov (400²) | near vacuum, very strong shock | 12/12 ok; shock radius within 0–2 % of the prediction, spread ≤ 1.7 % | – | PASS |
| Mach 80 jet (800 × 400) | positivity at high Mach | 11/12 ok; VK p = −9.5 | VK still p = −45 | PASS (MOOD), FAIL (VK) |
| Mach 2000 jet (640 × 320) | positivity at extreme Mach | MUSCL1 only; all K ≥ 2 end with p < 0 | all MOOD ok (min p 2e-4 … 3e-2); VK p = −2.7e4 | PASS for MOOD in Float64 only, FAIL for VK |
| Gresho (200², t = 3) | smooth vortex, accuracy | 12/12 ok, but MOOD loses 50–90 % of the vortex (L1(ρu) 0.08–0.15 of 0.16), VK 0.021 | – | PASS (runs), poor accuracy with MOOD |
| Gresho, Mach 0.01 (128², t = 1) | low-Mach dissipation | 12/12 ok, vortex gone for every method (L1(ρu) ≈ 0.13–0.16 of 0.16) | – | FAIL (accuracy), as expected without low-Mach fix |
| Kelvin–Helmholtz (512², t = 2) | numerical dissipation | 12/12 ok; mode grows only with VK (0.01 → 0.20), stays (EPD1) or decays (EPD2, MUSCL1) | – | FAIL with MOOD (instability suppressed) |
| Shock–vortex (800 × 400) | shock + smooth vortex | 12/12 ok | – | PASS |

## 1. Robustness: everything runs, positivity fails only in Float32 or with VK

- All 168 GPU runs reach t_end without NaN. The only failures are negative pressures at t_end:
  - **Float32 round-off** where the pressure is a tiny difference of large energies: the cold inflow of Noh (p0 = 1e-6
    against ½ρ|u|² = 0.5; min p −1.2e-7 … −8e-8, all MLS 3/4 runs and K = 4 EPD2 MLS 2) and the Mach 2000 jet (energy
    density 1.6e6, so Float32 resolves the pressure only to about ±0.1 there; min p −1.3 … −28, more than one rounding
    error, so it accumulates over the steps). **All of them are positive in
    Float64** (Noh: min p 4e-10 … 2e-8; Mach 2000: 2e-4 … 3e-2, table in §6). MOOD's admissibility check works; it
    cannot see a pressure that rounds to a negative number.
  - **The VK limiter is not positivity preserving**: Noh (p = −1.5e-2), Mach 80 (−9.5 / −45 in Float64) and Mach 2000
    (−3.3e4 / −2.7e4) fail in both precisions. VK has no a-posteriori safeguard; MOOD runs the same problems.
- MUSCL1 passes every test in Float32, including Mach 2000.
- **Library finding (HyperCloudCore 0.7.8): the CFL estimate ignores boundary states.** `_cfl_body!` skips boundary
  particles and uses each particle's own state, so a fast inflow (the jets) is not seen until it is inside: the first
  step of the jets was ~25× (Mach 80) / ~700× (Mach 2000) too large and produced p ≈ −60 / −6e7 after one step. The jet
  configs work around it by starting the jet in a thin strip x < 0.005 at the inlet. Proposed fix: take the maximum
  wave speed over the neighbourhood (including boundary particles) in the CFL estimate.
- Cost: frozen MLS 3/4 needs 3× the steps (the CFL follows the degree-3/4 divergence) and MOOD runs are 3–10× slower per
  step than MUSCL1/VK on the GPU (up to 9.7 h for the implosion with K = 4, MLS 4, EPD2).

## 2. Accuracy on smooth flows: MOOD is very dissipative on shear and vortex flows

- **Kelvin–Helmholtz** (McNally et al. 2012 smooth setup, seeded mode amplitude M(0) = 0.01; physically it grows to
  O(0.1) by t ≈ 1.5): M(2) = 0.196 with VK, 0.0055 … 0.011 with EPD1, 1.4e-4 … 2.4e-3 with EPD2, 4.1e-4 with MUSCL1.
  MOOD suppresses the instability about as much as first order.
- **Gresho vortex** (exact steady solution; the L1 norm of ρu of the exact solution is 0.16): L1(ρu) at t = 3 is 0.021
  with VK, 0.08–0.11 with EPD1, 0.13–0.15 with EPD2 and 0.157 with MUSCL1 (the vortex is dissipated almost completely).
  The velocity profile has kinks at r = 0.2 and 0.4, where the DMP triggers legitimately, but EPD2 at K ≥ 3 is no better
  than K = 2.
- **Low Mach (Ma = 0.01)**: every method, also VK, ends at L1(ρu) ≈ 0.13–0.16, i.e. the Rusanov-type dissipation (∝ the
  sound speed) destroys the vortex within one rotation. Expected without a low-Mach correction (the PDF's point).

Diagnosis (`tools/mood_activity.jl`, Kelvin–Helmholtz at 256², t ≤ 0.5): MOOD drops **2–3 % of the particles to order 1
in every step with EPD2 at K = 2, and 7–9 % below K (5–6 % at order 1) with K = 4 / MLS 4**, in Float32 and in
Float64 alike (not round-off). The dropped particles lie in the smoothed shear layers (median distance to the interface
0.03–0.05, the layer width is L = 0.025), not at the extrema of the seed. Changing the U2 parameters does not change
this (curvature ratio 0, plateau 1e-2, δ = 1e-2 Δx²: still 6–7 % below K). My guess (not verified): the curvature
*sign* test of U2 rejects these particles, because the McNally profile changes the sign of its curvature at the
interface (inflection), so every DMP violation near the layer is classified as non-smooth. Next step: log which test
rejects (DMP, plateau, sign, ratio) per particle in HyperCloudCore.

## 3. Symmetry (implosion)

On a lattice the x↔y mirror image of every particle is a particle, so the symmetry error Σ|ρ(x, y) − ρ(y, x)| / Σρ is
exact: 4.1e-7 (MUSCL1), 3.0e-6 (VK), 1.3e-4 … 8.8e-4 (MOOD, growing with K and with MLS = K). The jet along the
diagonal forms in all runs (min ρ 0.40–0.48 for K ≥ 2 vs 0.70 for MUSCL1). MOOD's order decisions are not symmetric
(per-particle decisions on neighbourhoods summed in different orders), which seeds the asymmetry.

## 4. Strong shocks with exact solutions

- **Sedov**: shock radius 0.979 (MUSCL1) … 1.000 (K = 4, MLS 4) against the predicted 1.001 from the deposited energy;
  angular spread 0.6 % (MUSCL1, VK) to 1.7 %; peak density 3.2 (MUSCL1) … 4.4 (exact 6, resolution limited at 400²).
  min ρ at the centre 0.006–0.011 for K ≥ 2 (0.043 with MUSCL1; near vacuum), all positive.
- **Noh**: L1(ρ) against the exact solution 0.67 (MUSCL1), 0.19–0.21 with EPD1, 0.21–0.36 with EPD2 (EPD2 is worse
  here), VK 0.23 (but negative pressure). The error is dominated by the shock and the centre (wall heating; not
  quantified separately).

## 5. Tables (Float32, GPU)

Generated with `julia --project=env_headless tools/bench_table.jl <configs>`.

### Riemann2D_C4
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.25 / 0.25 | 0 | 0.5065 | 3.500e-01 | 3239 | 107 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5057 | 3.482e-01 | 3286 | 124 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5052 | 3.487e-01 | 3279 | 301 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5060 | 3.495e-01 | 3257 | 867 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5055 | 3.476e-01 | 3282 | 1196 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5061 | 3.497e-01 | 3261 | 945 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.25 / 0.25 | 0 | 0.5063 | 3.410e-01 | 10137 | 5399 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.25 / 0.25 | 0 | 0.5062 | 3.497e-01 | 10139 | 6291 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5063 | 3.498e-01 | 3264 | 1328 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5064 | 3.500e-01 | 3245 | 3544 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.25 / 0.25 | 0 | 0.5058 | 3.493e-01 | 10133 | 10535 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.25 / 0.25 | 0 | 0.5043 | 3.468e-01 | 10146 | 5043 |  |
### Riemann2D_C6
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.3 / 0.3 | 0 | 0.4870 | 2.621e-01 | 3401 | 121 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.3 / 0.3 | 0 | 0.2583 | 1.254e-01 | 3399 | 138 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.3 / 0.3 | 0 | 0.2410 | 1.236e-01 | 3387 | 533 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.3 / 0.3 | 0 | 0.2903 | 1.469e-01 | 3387 | 786 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.3 / 0.3 | 0 | 0.2408 | 1.143e-01 | 3389 | 1082 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.3 / 0.3 | 0 | 0.2943 | 1.425e-01 | 3388 | 666 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.3 / 0.3 | 0 | 0.2377 | 1.205e-01 | 10343 | 5120 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.3 / 0.3 | 0 | 0.2552 | 1.275e-01 | 10336 | 1951 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.3 / 0.3 | 0 | 0.2584 | 1.334e-01 | 3389 | 2739 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.3 / 0.3 | 0 | 0.3040 | 1.543e-01 | 3393 | 3412 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.3 / 0.3 | 0 | 0.2632 | 1.376e-01 | 10340 | 4548 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.3 / 0.3 | 0 | 0.2684 | 1.338e-01 | 10339 | 13256 |  |
### Riemann2D_C12
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.25 / 0.25 | 0 | 0.5313 | 4.000e-01 | 2928 | 108 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5311 | 3.997e-01 | 3274 | 134 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5308 | 3.995e-01 | 3208 | 556 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5307 | 3.994e-01 | 3073 | 735 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5311 | 3.998e-01 | 3220 | 1094 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5310 | 3.997e-01 | 3077 | 686 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.25 / 0.25 | 0 | 0.5311 | 3.998e-01 | 9150 | 1746 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.25 / 0.25 | 0 | 0.5312 | 3.999e-01 | 8983 | 5416 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5312 | 3.999e-01 | 3242 | 1293 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.25 / 0.25 | 0 | 0.5313 | 4.000e-01 | 3084 | 1261 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.25 / 0.25 | 0 | 0.5311 | 3.998e-01 | 9169 | 3930 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.25 / 0.25 | 0 | 0.5308 | 3.995e-01 | 8996 | 10945 |  |
### Riemann2D_C15
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.2 / 0.2 | 0 | 0.5120 | 3.993e-01 | 1869 | 76 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5110 | 3.406e-01 | 1892 | 86 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5113 | 3.463e-01 | 1897 | 359 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5111 | 3.850e-01 | 1879 | 466 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5097 | 3.471e-01 | 1892 | 677 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5112 | 3.838e-01 | 1873 | 383 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.2 / 0.2 | 0 | 0.5108 | 3.482e-01 | 5879 | 1186 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.2 / 0.2 | 0 | 0.5083 | 3.701e-01 | 5878 | 1210 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5122 | 3.450e-01 | 1885 | 770 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.2 / 0.2 | 0 | 0.5118 | 3.867e-01 | 1872 | 1908 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.2 / 0.2 | 0 | 0.5111 | 3.548e-01 | 5878 | 2467 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.2 / 0.2 | 0 | 0.5074 | 3.713e-01 | 5880 | 3176 |  |
### Riemann2D_C17
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.3 / 0.3 | 0 | 0.5108 | 3.555e-01 | 3152 | 114 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.3 / 0.3 | 0 | 0.3260 | 1.410e-01 | 3313 | 135 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.3 / 0.3 | 0 | 0.3286 | 1.759e-01 | 3343 | 585 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.3 / 0.3 | 0 | 0.4412 | 2.411e-01 | 3291 | 801 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.3 / 0.3 | 0 | 0.3431 | 1.797e-01 | 3335 | 1161 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.3 / 0.3 | 0 | 0.4557 | 2.431e-01 | 3289 | 617 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.3 / 0.3 | 0 | 0.3647 | 1.934e-01 | 9699 | 1872 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.3 / 0.3 | 0 | 0.4140 | 2.237e-01 | 9666 | 5968 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.3 / 0.3 | 0 | 0.3580 | 1.853e-01 | 3339 | 2726 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.3 / 0.3 | 0 | 0.4570 | 2.484e-01 | 3267 | 3214 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.3 / 0.3 | 0 | 0.3854 | 2.018e-01 | 9669 | 4245 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.3 / 0.3 | 0 | 0.4187 | 2.218e-01 | 9648 | 4395 |  |
### Implosion
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | symmetry error (final) |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 2.5 / 2.5 | 0 | 0.7000 | 7.859e-01 | 32384 | 304 | 4.058e-07 |
| RK2MUSCL2(VK,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4102 | 7.499e-01 | 39443 | 414 | 2.993e-06 |
| RK2MUSCL2(EPD1,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4112 | 7.357e-01 | 40304 | 2365 | 1.322e-04 |
| RK2MUSCL2(EPD2,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4660 | 7.545e-01 | 39051 | 2912 | 3.337e-04 |
| RK3MUSCL3(EPD1,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4292 | 7.503e-01 | 39889 | 4782 | 2.284e-04 |
| RK3MUSCL3(EPD2,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4772 | 7.597e-01 | 38777 | 2628 | 4.552e-04 |
| RK3MUSCL3(EPD1,AR,MLS3) | 2.5 / 2.5 | 0 | 0.3972 | 7.489e-01 | 105957 | 15700 | 3.629e-04 |
| RK3MUSCL3(EPD2,AR,MLS3) | 2.5 / 2.5 | 0 | 0.3986 | 7.485e-01 | 105681 | 17919 | 5.746e-04 |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4156 | 7.533e-01 | 40021 | 4661 | 5.811e-04 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 2.5 / 2.5 | 0 | 0.4459 | 7.585e-01 | 39098 | 12086 | 5.279e-04 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 2.5 / 2.5 | 0 | 0.3979 | 7.548e-01 | 105524 | 13253 | 7.970e-04 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 2.5 / 2.5 | 0 | 0.4013 | 7.539e-01 | 105536 | 35102 | 8.835e-04 |
### Noh
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | L1(ρ) / L2(ρ) / L1(ρu) vs exact |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.6 / 0.6 | 0 | 1.4243 | 1.240e-05 | 818 | 27 | 6.669e-01 / 1.647e+00 / 9.263e-02 |
| RK2MUSCL2(VK,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4257 | -1.517e-02 | 802 | 28 | 2.298e-01 / 7.849e-01 / 5.535e-02 |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4258 | 3.576e-07 | 789 | 67 | 2.128e-01 / 7.911e-01 / 3.620e-02 |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4258 | 3.179e-07 | 800 | 80 | 3.476e-01 / 1.076e+00 / 5.190e-02 |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4260 | 3.974e-08 | 781 | 106 | 2.113e-01 / 7.774e-01 / 3.681e-02 |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4260 | 7.947e-08 | 800 | 56 | 3.491e-01 / 1.070e+00 / 5.163e-02 |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.6 / 0.6 | 0 | 1.4253 | -1.192e-07 | 2254 | 328 | 1.976e-01 / 7.366e-01 / 3.674e-02 |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.6 / 0.6 | 0 | 1.4253 | -7.947e-08 | 2264 | 393 | 2.119e-01 / 7.704e-01 / 3.835e-02 |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4259 | 3.974e-08 | 795 | 183 | 1.907e-01 / 7.025e-01 / 4.378e-02 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4259 | 0.000e+00 | 805 | 239 | 3.637e-01 / 1.053e+00 / 6.753e-02 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.6 / 0.6 | 0 | 1.4252 | -7.947e-08 | 2260 | 305 | 1.950e-01 / 6.859e-01 / 4.979e-02 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.6 / 0.6 | 0 | 1.4252 | -7.947e-08 | 2256 | 300 | 2.310e-01 / 7.555e-01 / 5.986e-02 |
### Sedov
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | R_shock / R_pred / spread / ρ_max |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 1 / 1 | 0 | 0.0432 | 1.000e-05 | 2716 | 46 | 0.979 / 1.001 / 0.006 / 3.23 |
| RK2MUSCL2(VK,AR,MLS2) | 1 / 1 | 0 | 0.0056 | 1.000e-05 | 6456 | 87 | 0.986 / 1.001 / 0.006 / 4.38 |
| RK2MUSCL2(EPD1,AR,MLS2) | 1 / 1 | 0 | 0.0064 | 6.506e-06 | 6030 | 339 | 0.985 / 1.001 / 0.011 / 4.14 |
| RK2MUSCL2(EPD2,AR,MLS2) | 1 / 1 | 0 | 0.0106 | 9.944e-06 | 4742 | 359 | 0.985 / 1.001 / 0.014 / 3.77 |
| RK3MUSCL3(EPD1,AR,MLS2) | 1 / 1 | 0 | 0.0069 | 3.382e-06 | 5841 | 594 | 0.985 / 1.001 / 0.010 / 4.17 |
| RK3MUSCL3(EPD2,AR,MLS2) | 1 / 1 | 0 | 0.0108 | 9.907e-06 | 4709 | 442 | 0.987 / 1.001 / 0.014 / 3.79 |
| RK3MUSCL3(EPD1,AR,MLS3) | 1 / 1 | 0 | 0.0062 | 6.150e-07 | 17829 | 2922 | 0.993 / 1.001 / 0.017 / 4.07 |
| RK3MUSCL3(EPD2,AR,MLS3) | 1 / 1 | 0 | 0.0066 | 4.408e-06 | 17381 | 3416 | 0.998 / 1.001 / 0.015 / 4.04 |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 1 / 1 | 0 | 0.0062 | 9.726e-06 | 6212 | 1371 | 0.986 / 1.001 / 0.010 / 4.28 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 1 / 1 | 0 | 0.0099 | 9.928e-06 | 4929 | 1292 | 0.987 / 1.001 / 0.016 / 3.82 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 1 / 1 | 0 | 0.0069 | 2.358e-07 | 17597 | 2520 | 0.999 / 1.001 / 0.017 / 4.30 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 1 / 1 | 0 | 0.0073 | 2.491e-07 | 17158 | 2543 | 1.000 / 1.001 / 0.013 / 4.15 |
### MachJet80
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.07 / 0.07 | 0 | 0.1424 | 4.123e-01 | 4958 | 99 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.07 / 0.07 | 0 | 0.1158 | -9.471e+00 | 5528 | 121 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.07 / 0.07 | 0 | 0.1197 | 3.255e-04 | 6308 | 630 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.07 / 0.07 | 0 | 0.1196 | 1.361e-01 | 6015 | 783 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.07 / 0.07 | 0 | 0.1159 | 6.510e-04 | 5684 | 577 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.07 / 0.07 | 0 | 0.1017 | 2.061e-01 | 5613 | 600 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.07 / 0.07 | 0 | 0.0136 | 5.315e-04 | 16882 | 2578 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.07 / 0.07 | 0 | 0.0488 | 1.790e-03 | 16799 | 5944 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.07 / 0.07 | 0 | 0.0977 | 1.302e-03 | 5530 | 2579 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.07 / 0.07 | 0 | 0.0999 | 2.020e-01 | 5475 | 2754 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.07 / 0.07 | 0 | 0.0294 | 6.510e-04 | 16774 | 4041 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.07 / 0.07 | 0 | 0.0078 | 1.059e-03 | 16322 | 11653 |  |
### MachJet2000
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.001 / 0.001 | 0 | 0.1716 | 2.500e-01 | 3165 | 58 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1357 | -3.330e+04 | 3526 | 67 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1181 | -1.283e+01 | 3883 | 322 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1415 | -8.083e+00 | 3713 | 382 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1191 | -9.750e+00 | 3597 | 337 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1331 | -1.392e+01 | 3559 | 717 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.001 / 0.001 | 0 | 0.0201 | -2.233e+01 | 10541 | 2323 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.001 / 0.001 | 0 | 0.0032 | -1.250e+00 | 10486 | 1226 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1137 | -2.833e+01 | 3516 | 1197 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1219 | -1.225e+01 | 3487 | 527 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.001 / 0.001 | 0 | 0.0022 | -1.200e+01 | 10894 | 4205 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.001 / 0.001 | 0 | 0.0383 | -1.250e+01 | 10348 | 5143 |  |
### Gresho
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | L1(ρ) / L2(ρ) / L1(ρu) vs exact |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 3 / 3 | 0 | 0.9973 | 5.721e+00 | 13601 | 62 | 1.161e-03 / 1.420e-03 / 1.570e-01 |
| RK2MUSCL2(VK,AR,MLS2) | 3 / 3 | 0 | 0.9907 | 5.022e+00 | 15170 | 69 | 1.108e-03 / 1.407e-03 / 2.138e-02 |
| RK2MUSCL2(EPD1,AR,MLS2) | 3 / 3 | 0 | 0.9940 | 5.370e+00 | 14781 | 246 | 4.383e-03 / 5.755e-03 / 8.101e-02 |
| RK2MUSCL2(EPD2,AR,MLS2) | 3 / 3 | 0 | 0.9931 | 5.703e+00 | 13843 | 144 | 4.055e-03 / 5.965e-03 / 1.438e-01 |
| RK3MUSCL3(EPD1,AR,MLS2) | 3 / 3 | 0 | 0.9943 | 5.420e+00 | 14690 | 232 | 4.550e-03 / 6.014e-03 / 8.802e-02 |
| RK3MUSCL3(EPD2,AR,MLS2) | 3 / 3 | 0 | 0.9937 | 5.707e+00 | 13798 | 720 | 3.795e-03 / 5.566e-03 / 1.459e-01 |
| RK3MUSCL3(EPD1,AR,MLS3) | 3 / 3 | 0 | 0.9944 | 5.442e+00 | 42191 | 949 | 4.506e-03 / 5.945e-03 / 9.315e-02 |
| RK3MUSCL3(EPD2,AR,MLS3) | 3 / 3 | 0 | 0.9928 | 5.659e+00 | 41292 | 2358 | 4.662e-03 / 6.316e-03 / 1.286e-01 |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 3 / 3 | 0 | 0.9942 | 5.493e+00 | 14576 | 570 | 4.659e-03 / 5.980e-03 / 9.865e-02 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 3 / 3 | 0 | 0.9937 | 5.709e+00 | 13779 | 1552 | 3.626e-03 / 5.275e-03 / 1.470e-01 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 3 / 3 | 0 | 0.9935 | 5.566e+00 | 41805 | 4604 | 4.871e-03 / 6.289e-03 / 1.100e-01 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 3 / 3 | 0 | 0.9929 | 5.681e+00 | 41162 | 5072 | 4.704e-03 / 6.750e-03 / 1.347e-01 |
### GreshoLowMach
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | L1(ρ) / L2(ρ) / L1(ρu) vs exact |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 97985 | 99 | 1.111e-06 / 1.523e-06 / 1.601e-01 |
| RK2MUSCL2(VK,AR,MLS2) | 1 / 1 | 0 | 0.9998 | 7.144e+03 | 98157 | 114 | 4.056e-05 / 6.572e-05 / 1.338e-01 |
| RK2MUSCL2(EPD1,AR,MLS2) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 98028 | 837 | 1.445e-05 / 1.476e-05 / 1.600e-01 |
| RK2MUSCL2(EPD2,AR,MLS2) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 97985 | 1344 | 6.090e-06 / 6.912e-06 / 1.601e-01 |
| RK3MUSCL3(EPD1,AR,MLS2) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 98027 | 1823 | 1.327e-05 / 1.396e-05 / 1.599e-01 |
| RK3MUSCL3(EPD2,AR,MLS2) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 97985 | 2726 | 6.035e-06 / 6.665e-06 / 1.600e-01 |
| RK3MUSCL3(EPD1,AR,MLS3) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 284386 | 6689 | 1.801e-05 / 1.834e-05 / 1.602e-01 |
| RK3MUSCL3(EPD2,AR,MLS3) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 284379 | 9360 | 1.321e-05 / 1.377e-05 / 1.604e-01 |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 98025 | 4418 | 1.396e-05 / 1.598e-05 / 1.600e-01 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 97985 | 3658 | 8.753e-06 / 1.023e-05 / 1.601e-01 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 1 / 1 | 0 | 0.9999 | 7.144e+03 | 284383 | 7395 | 2.206e-05 / 2.817e-05 / 1.603e-01 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 1 / 1 | 0 | 1.0000 | 7.144e+03 | 284378 | 18956 | 1.630e-05 / 1.915e-05 / 1.606e-01 |
### KelvinHelmholtz
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | KH mode amplitude M(t_end) (M(0) = 0.01) |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 2 / 2 | 0 | 1.0823 | 2.550e+00 | 17856 | 480 | 4.093e-04 |
| RK2MUSCL2(VK,AR,MLS2) | 2 / 2 | 0 | 0.9944 | 2.203e+00 | 18313 | 556 | 1.964e-01 |
| RK2MUSCL2(EPD1,AR,MLS2) | 2 / 2 | 0 | 1.0009 | 2.480e+00 | 18197 | 514 | 1.075e-02 |
| RK2MUSCL2(EPD2,AR,MLS2) | 2 / 2 | 0 | 1.0142 | 2.524e+00 | 18187 | 2363 | 1.360e-04 |
| RK3MUSCL3(EPD1,AR,MLS2) | 2 / 2 | 0 | 1.0009 | 2.478e+00 | 18196 | 2799 | 1.118e-02 |
| RK3MUSCL3(EPD2,AR,MLS2) | 2 / 2 | 0 | 1.0155 | 2.524e+00 | 18166 | 1523 | 2.732e-04 |
| RK3MUSCL3(EPD1,AR,MLS3) | 2 / 2 | 0 | 1.0026 | 2.486e+00 | 57238 | 12353 | 8.490e-03 |
| RK3MUSCL3(EPD2,AR,MLS3) | 2 / 2 | 0 | 1.0061 | 2.504e+00 | 57232 | 6116 | 2.156e-03 |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 2 / 2 | 0 | 1.0032 | 2.497e+00 | 18200 | 8000 | 8.538e-03 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 2 / 2 | 0 | 1.0222 | 2.527e+00 | 18154 | 3660 | 6.201e-04 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 2 / 2 | 0 | 1.0045 | 2.501e+00 | 57241 | 15568 | 5.510e-03 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 2 / 2 | 0 | 1.0073 | 2.508e+00 | 57230 | 13042 | 2.379e-03 |
### ShockVortex
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1 | 0.6 / 0.6 | 0 | 0.9948 | 9.927e-01 | 3548 | 76 |  |
| RK2MUSCL2(VK,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9909 | 9.873e-01 | 3689 | 87 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9901 | 9.862e-01 | 3690 | 324 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9901 | 9.862e-01 | 3675 | 462 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9898 | 9.857e-01 | 3694 | 625 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9898 | 9.857e-01 | 3681 | 362 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.6 / 0.6 | 0 | 0.9892 | 9.849e-01 | 10974 | 1003 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.6 / 0.6 | 0 | 0.9892 | 9.849e-01 | 10956 | 1081 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9894 | 9.852e-01 | 3693 | 719 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 0.9895 | 9.853e-01 | 3683 | 1814 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.6 / 0.6 | 0 | 0.9889 | 9.845e-01 | 10974 | 5663 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.6 / 0.6 | 0 | 0.9887 | 9.843e-01 | 10951 | 6630 |  |

## 6. Float64 reruns (CPU, local)

### Noh_Float64
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | L1(ρ) / L2(ρ) / L1(ρu) vs exact |
|---|---|---|---|---|---|---|---|
| RK2MUSCL2(VK,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4258 | -1.567e-02 | 801 | 17 | 2.301e-01 / 7.886e-01 / 5.557e-02 |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.6 / 0.6 | 0 | 1.4253 | 1.037e-09 | 2203 | 135 | 2.028e-01 / 7.498e-01 / 3.843e-02 |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.6 / 0.6 | 0 | 1.4253 | 3.052e-09 | 2225 | 163 | 2.178e-01 / 7.835e-01 / 4.143e-02 |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.6 / 0.6 | 0 | 1.4260 | 2.134e-08 | 802 | 106 | 3.712e-01 / 1.065e+00 / 7.251e-02 |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.6 / 0.6 | 0 | 1.4252 | 3.576e-10 | 2219 | 302 | 2.000e-01 / 7.071e-01 / 5.113e-02 |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.6 / 0.6 | 0 | 1.4252 | 6.117e-10 | 2221 | 339 | 2.371e-01 / 7.784e-01 / 6.019e-02 |
### MachJet80_Float64
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL2(VK,AR,MLS2) | 0.07 / 0.07 | 0 | 0.1659 | -4.512e+01 | 5527 | 261 |  |
### MachJet2000_Float64
| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] |  |
|---|---|---|---|---|---|---|---|
| RK2MUSCL2(VK,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1672 | -2.703e+04 | 3514 | 107 |  |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1574 | 6.840e-04 | 3866 | 164 |  |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1658 | 1.930e-02 | 3706 | 190 |  |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1157 | 3.320e-04 | 3584 | 289 |  |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1583 | 2.651e-02 | 3550 | 338 |  |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.001 / 0.001 | 0 | 0.0239 | 4.954e-04 | 10636 | 1160 |  |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.001 / 0.001 | 0 | 0.0331 | 1.315e-03 | 10618 | 1393 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1079 | 4.834e-04 | 3505 | 698 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.001 / 0.001 | 0 | 0.1557 | 2.360e-04 | 3475 | 657 |  |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.001 / 0.001 | 0 | 0.0586 | 2.668e-04 | 10731 | 3439 |  |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.001 / 0.001 | 0 | 0.0373 | 1.544e-03 | 10431 | 2546 |  |

## Verdict

**INCONCLUSIVE overall**: the code is robust (no NaN, no blow-up in 168 runs; positivity holds for MOOD in Float64 even
at Mach 2000), but in single precision the cold-gas and hypersonic cases lose positivity through round-off, the VK
limiter is not positivity preserving, and MOOD U2 is too dissipative on smooth shear and vortex flows (Kelvin–Helmholtz,
Gresho). Proposed library work: (1) CFL estimate including neighbour/boundary states, (2) find which U2 test rejects
smooth shear layers, (3) a positivity limiter for VK (or VK + MOOD admissibility), (4) a low-Mach correction of the
dissipation.
