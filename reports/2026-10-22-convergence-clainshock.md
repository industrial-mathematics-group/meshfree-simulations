# 2026-10-22: Convergence sweeps (new package stack) and ClainShock MLS-order runs

Branch `claude/convergence` (from `main` 728fd1c), cluster jobs through the gateway. Package stack of the runs
(re-resolved Manifest): **HyperCloud 0.5.0, HyperCloudCore 0.7.0, PDEStudioCore 0.4.1** (the `__all` sweep
suffix), CUDACore from CUDA.jl `main`. The ClainShock runs are therefore not bit-comparable with the
Mach3Step/DMR runs of REPORT_2026-10-21.md (HyperCloudCore 0.6.2).

## 1. Convergence (Float64, CFL 0.2, R = 5.5, RK4, analytical reference; regular lattice, jitter 0 in all three configs)

Jobs: Convergence1D 24544690–693 (CPU, 32 cores, ~1.5 min each); 2-D sweeps N ≤ 200 on CPU
(24544702–705, 24544823/831/844/857, 1–3 min each) and N ≥ 251 per method on a V100 (ConvergenceAdvection
24544694–697, ConvergenceBurgers 24544709–711/716). All 84 + 2 × 64 runs completed; mass drift 1e-16 everywhere.
EOC = log(e_{i-1}/e_i) / log(N_i/N_{i-1}) of the L1 error of u₁, last five refinements.

| Experiment | MUSCL2 | MUSCL3 | MUSCL4 | MUSCL5 |
|---|---|---|---|---|
| 1-D advection, N = 32 … 3163, t = 10 | 2.00 (1.0e-3) | **3.00** (1.6e-5) | 4.00 (2.1e-8) | 5.00 (5.3e-10) |
| 2-D advection, N² ≤ 1000², t = 1 | 2.06 (2.7e-3) | 3.00 (1.9e-4) | **4.1, then floor 5.3e-6 from N = 631** | 4.98 up to N = 252, **then floor 2.1e-5 from N = 399** |
| 2-D Burgers, N² ≤ 1000², t = 0.5 | 2.12 (8.2e-4) | 2.94 (8.8e-5) | 4.2 → 3.2 at N = 1000 (2.0e-6) | 4.3 at N ≈ 400, **then floor 3.7e-6** (EOC 2.3, 0.6) |

(in brackets: L1 error at the finest grid)

- **MUSCL3 is third order now** in 1-D and 2-D (2.0 in the 2026-10-17 study): the degree-3 divergence fit of
  HyperCloudCore ≥ 0.6.2 for odd orders.
- **MUSCL4 on 2-D advection stops converging at N ≈ 600 on the GPU**: errors 9.9e-6 (502²), 5.40e-6 (631²),
  5.28e-6 (795²), 5.49e-6 (1000²). The 2026-10-17 study reached 6e-7 at 1000² with order 4.03. Burgers shows
  the start of the same floor (EOC 4.2 → 3.2 at 1000²). **Cause: the CUDA path, not the numerics**, see §1b.
- MUSCL5 reaches 5.0 in 1-D. In 2-D it converges at 4.6–5.0 up to N ≈ 250 (advection: 7.5e-5 at 252²) and
  then stops dead: 3.0e-5, 2.2e-5, 2.1e-5, 2.1e-5, 2.1e-5, 2.1e-5 (317² … 1000²). Burgers: peak 4.3 at
  N = 317–502, then 7.3e-6 (631²), 4.3e-6 (795²), 3.7e-6 (1000²).
- So every 2-D sweep of order ≥ 4 that ran on the V100 (N ≥ 251) ends on a flat floor: MUSCL4 5e-6 / ~2e-6,
  MUSCL5 2.1e-5 / 3.7e-6 (advection / Burgers); the 1-D sweeps and the N ≤ 200 runs (all CPU) do not.

### 1b. The floor is a GPU artefact (CPU reference runs, same configuration)

The N ≥ 251 advection sweeps repeated on CPU nodes (32 cores, Float64, identical parameters; jobs 24559671
MUSCL4, 24559672/673 MUSCL5), next to the V100 values (24544696/697, reproduced bit-identically by 24558179/180):

| N | MUSCL4 CPU | MUSCL4 V100 | MUSCL5 CPU | MUSCL5 V100 |
|---|---|---|---|---|
| 252 | 1.638e-4 | 1.634e-4 | 7.313e-5 | 7.503e-5 |
| 317 | 6.298e-5 | 6.261e-5 | 2.326e-5 | 3.016e-5 |
| 399 | 2.448e-5 | 2.429e-5 | 7.374e-6 | 2.232e-5 |
| 502 | 9.615e-6 | 9.917e-6 | 2.341e-6 | 2.117e-5 |
| 631 | 3.813e-6 | 5.399e-6 | 7.464e-7 | 2.103e-5 |
| 795 | 1.503e-6 | 5.283e-6 | 2.352e-7 | 2.102e-5 |
| 1000 | **5.979e-7** | 5.489e-6 | **7.470e-8** | 2.102e-5 |

On the CPU both schemes keep their order (MUSCL4: 3.9–4.0 up to 1000², 6.0e-7 as in the 2026-10-17 study;
MUSCL5: 5.0 all the way to 1000², 7.5e-8). On the V100 the errors leave the CPU values already at N = 252 (MUSCL5: 7.50e-5 vs
7.31e-5) and freeze at 2.1e-5 / 5e-6. Divided by the box area (100) the floor is a mean error of ~2e-7, i.e.
Float32 epsilon, although T = Float64 (checked in the stored states and positions). **Guess:** a Float32
operation or constant in the CUDA kernels/extension (weight function, divergence weights, or the neighbour
slot data), not in HyperCloudCore's interpolation (which the CPU runs exercise unchanged). Float32 runs of
the same N = 399 MUSCL5 case (jobs 24561401 CPU, 24561402 V100): CPU 3.95e-5, GPU 3.43e-5, i.e. the GPU's
*Float64* result (2.23e-5) sits at the Float32 error level, while the CPU Float64 result is 7.37e-6. So
part of the device computation runs in single precision although T = Float64 (not all of it: the GPU
Float64 value is still below the Float32 ones).

The GPU-based ClainShock/Mach3Step/DMR runs are Float32 anyway, so they are not affected by this beyond their
own precision; the convergence conclusions for MUSCL4/5 in 2-D must be taken from the CPU runs.

## 2. ClainShock 800² (Lax–Liu 3), CFL 0.8, Float32, t = 0.3, MLS-order variants

Jobs 24557467–472, 24557550, 24557635, 24557812 (V100, dgx001 / gpu004). "mass" = final mass relative to the
start (the outer states flow in through the boundaries), Δmass relative to the final mass of RK2MUSCL2(EPD1).

| Method | MLS | min p (at) | mass | Δmass | steps | runtime |
|---|---|---|---|---|---|---|
| RK2MUSCL2(EPD1) (reference) | 2 | 0.018 (0.40, 0.24) | 0.379 | – | 1808 | 2.6 min |
| RK3MUSCL3(EPD1) | 0 | **0.0028** (0.22, 0.38) | **0.424** | **+3.3 %** | 5445 | 18 min |
| RK3MUSCL3(EPD2) | 0 | 0.010 | 0.372 | −0.5 % | 5393 | 22 min |
| RK3MUSCL3(EPD1) | 2 | 0.015 | 0.377 | −0.1 % | 1812 | 6 min |
| RK3MUSCL3(EPD1) | **3** | **< 0 at t = 0.0019**, x = (0.002, 0.500) | – | – | 36 | – |
| RK4MUSCL4(EPD1) | 0 | 0.0051 | 0.395 | +1.2 % | 5589 | 37 min |
| RK4MUSCL4(EPD2) | 0 | 0.023 | 0.368 | −0.8 % | 5402 | 45 min |
| RK4MUSCL4(EPD1) | 2 | 0.023 | 0.372 | −0.5 % | 1825 | 11 min |
| RK4MUSCL4(EPD1) | **4** | **< 0 at t = 0.0018**, x = (0.001, 0.500) | – | – | 34 | – |

- **The fixed orders 3 and 4 fail here too, and there is no wall**: both blow up within 35 steps at the left
  domain boundary at y = 0.5, where the initial discontinuity between the two left quadrants meets the boundary.
  So the failure needs a discontinuity at a boundary (incomplete stencils, permanently lowered orders next to
  full-degree divergences), not a stagnation point. Same picture as Mach3Step (corner) and DMR (shock at the
  bottom wall).
- MLS 0 + EPD1 again gains mass (+3.3 % MUSCL3, +1.2 % MUSCL4) and has the lowest pressures; EPD2 and the fixed
  MLS 2 are within −0.1 … −0.8 % of the reference. The conservative variants spread more here (0.372–0.379)
  than on Mach3Step/DMR (1e-3); with inflow through all four boundaries the final mass also depends on how
  each scheme resolves the boundary states, so Δmass is a cruder measure on this problem.
- Cost at 800² on a V100: MUSCL2 2.6 min, MUSCL3 MLS 2 6 min, MUSCL3 MLS 0 18–22 min, MUSCL4 MLS 0 37–45 min.

## 2b. Float64 vs Float32 on CPU: the MLS-0 conservation error is not a precision artefact (2026-10-23)

Because the CUDA path is suspect in Float64 (§1b), the MUSCL3 Euler cases were rerun on the **CPU nodes**
(32 cores) in both precisions, CFL 0.8: ClainShock 800², Mach3Step 600×200, DMR 960×300; MUSCL2 EPD1 as the
reference, MUSCL3 with MLS 0 (EPD1, EPD2) and fixed MLS 2. Δmass relative to the final mass of the reference
in the same precision (Float32 and Float64 are different point clouds, the jitter is drawn in T). Jobs
24562745–755, 24562851 (Float64), 24562854–24564109 (Float32).

| Case | Method | Float64: min p / mass / Δ | Float32: min p / mass / Δ | Float32 on V100 (earlier) |
|---|---|---|---|---|
| ClainShock 800² | MUSCL2 EPD1 | 0.017 / 0.381 / – | 0.019 / 0.379 / – | 0.018 / 0.379 |
| | MUSCL3 EPD1, MLS 0 | 0.0025 / 0.419 / **+2.8 %** | 0.0034 / 0.422 / **+3.1 %** | 0.0028 / 0.424 |
| | MUSCL3 EPD2, MLS 0 | 0.0071 / 0.367 / −1.0 % | 0.020 / 0.372 / −0.5 % | 0.010 / 0.372 |
| | MUSCL3 EPD1, MLS 2 | 0.017 / 0.378 / −0.2 % | 0.015 / 0.377 / −0.1 % | 0.015 / 0.377 |
| Mach3Step 600×200 | MUSCL2 EPD1 | 0.445 / 0.795 / – | 0.334 / 0.773 / – | (CFL 0.4: 0.771) |
| | MUSCL3 EPD1, MLS 0 | 0.148 / 0.896 / **+5.6 %** | 0.196 / 0.922 / **+8.4 %** | 0.294 / 0.907 |
| | MUSCL3 EPD2, MLS 0 | 0.380 / 0.745 / −2.8 % | 0.179 / 0.746 / −1.5 % | 0.223 / 0.747 |
| | MUSCL3 EPD1, MLS 2 | 0.606 / 0.791 / −0.2 % | 0.335 / 0.770 / −0.2 % | 0.337 / 0.770 |
| DMR 960×300 | MUSCL2 EPD1 | 1.000 / 3.05 / – | 1.000 / 3.05 / – | 1.000 / 3.05 |
| | MUSCL3 EPD1, MLS 0 | 0.991 / 3.19 / **+3.5 %** | 0.990 / 3.19 / **+3.5 %** | 1.000 / 3.19 |
| | MUSCL3 EPD2, MLS 0 | 0.936 / 3.15 / +2.5 % | 0.985 / 3.15 / +2.5 % | – |
| | MUSCL3 EPD1, MLS 2 | 1.000 / 3.05 / 0 | 1.000 / 3.05 / 0 | – |

Higher resolutions, same setup (jobs 24562865–24564077):

| Case | Method | Float64: min p / mass / Δ | Float32: min p / mass / Δ |
|---|---|---|---|
| Mach3Step 900×300 | MUSCL2 EPD1 | 0.185 / 0.799 / – | 0.185 / 0.803 / – |
| | MUSCL3 EPD1, MLS 0 | 0.191 / 0.882 / **+4.6 %** | **0.00025** / 0.990 / **+10.4 %** |
| | MUSCL3 EPD2, MLS 0 | 0.115 / 0.735 / −3.6 % | 0.039 / 0.745 / −3.2 % |
| | MUSCL3 EPD1, MLS 2 | 0.193 / 0.798 / −0.1 % | 0.193 / 0.804 / +0.1 % |
| DMR 1440×450 | MUSCL2 EPD1 | 1.000 / 3.04 / – | 1.000 / 3.05 / – |
| | MUSCL3 EPD1, MLS 0 | 0.997 / 3.19 / +3.7 % | 0.963 / 3.19 / +3.5 % |
| | MUSCL3 EPD2, MLS 0 | 0.969 / 3.15 / +2.7 % | 0.966 / 3.15 / +2.5 % |
| | MUSCL3 EPD1, MLS 2 | 1.000 / 3.05 / +0.2 % | 1.000 / 3.05 / 0 |

- **Precision does not change the picture**: MLS 0 + EPD1 gains mass in Float64 too (ClainShock +2.8 %,
  DMR +3.5/+3.7 %, Mach3Step +5.6/+4.6 %), EPD2 is within ±4 %, the fixed MLS 2 within 0.2 % of the
  reference, on all three problems and both resolutions. DMR is identical to three digits in both
  precisions. On Mach3Step, Float64 does make MLS 0 + EPD1 *better behaved* (900×300: +4.6 % instead of
  +10.4 %, min p 0.19 instead of 2.5e-4, which was close to failing), so part of its erratic behaviour on
  that problem is Float32 sensitivity at the step corner, but the conservation error itself remains.
- Float32 on CPU vs Float32 on the V100 (same point cloud): identical for the conservative variants and
  for EPD2 (mass to 3 digits, step counts ±1), but **not for MLS 0 + EPD1** (Mach3Step: mass 0.922 vs
  0.907, min p 0.196 vs 0.294; ClainShock 0.422 vs 0.424). The single-precision effect in the device
  path (§1b) therefore also shows in Float32 runs, through the MOOD decisions of the most sensitive
  variant; the fixed-order results are reproducible across CPU and GPU.
- Cost on 32 EPYC cores: MUSCL2 3–13 min, MUSCL3 MLS 2 5–27 min, MUSCL3 MLS 0 23–140 min (Float64 ≈
  1.2–1.5 × Float32); Mach3Step 900×300 and DMR 1440×450 are the largest sizes that fit the 300 min limit
  for MLS 0 in Float64.

## 3. Verdict

- Convergence: PASS for all orders on CPU (2.0 / 3.0 / 4.0 / 5.0); **FAIL for the GPU path at orders 4–5 in
  2-D** (error floor ~2e-7 mean, Float32-like, while T = Float64). Needs a fix in the CUDA path.
- ClainShock: fixed MLS order 3/4 FAIL (boundary discontinuity); all MLS 0 / MLS 2 variants PASS (positive
  pressure, t = 0.3), with the known conservation error of MLS 0 + EPD1.
- Float64 CPU reruns (§2b): the MLS-0 conservation error is real, not a Float32 or GPU artefact; PASS for
  the fixed-order variants in both precisions.
