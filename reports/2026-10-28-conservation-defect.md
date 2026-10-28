# 2026-10-28: Conservation defect through the PDEStudioCore pipeline (paper section "conservation defect")

Recreation of the conservation-defect experiments (floor of the shock test, frozen/coupled divergence degree,
irregularity, support, weight, smooth 1-D/2-D) as `SimulationConfig`s, instead of the 2026-10-26 job scripts
(`cluster/run_mass.sh` with `SIM_SET`/`SIM_VARY` overrides) and the post-hoc `tools/mass_table.jl`.

- Branch `claude/defect-configs`. HyperCloud 0.5.1, HyperCloudCore 0.7.5, PDEStudioCore 0.4.2, Float64, CPU.
- **2374 runs, all local** (12 threads, tasks in parallel, each run serial; bitwise identical to a serial run
  on a smoke subset): floor 3.2 min, shock 18.8, jitter 2.7, range 2.4, α 1.5, smooth 1-D 3.0, smooth 2-D 21.9 min.
- Data: `results/defect/` (gitignored). Full tables: `reports/2026-10-28-conservation-defect-tables.md`.
- Rerun or retabulate:
  `julia -t 12 --project=env_headless tools/run_defect.jl defect_floor defect_shock …`
  (`DEFECT_SKIP_RUN=1` only tabulates). Each config also works with `@sim defect_shock`.

## Setup

`experiments/ConservationDefect/`: one config per paragraph (`defect_floor`, `defect_shock`, `defect_jitter`,
`defect_range`, `defect_alpha`, `defect_smooth1d`, `defect_smooth2d`) plus a shared `defect_setup.jl`.
- Shock problem and methods are unchanged: `conservation_burgers` and `methods_review.jl`.
- The sweep is `Grid_Ns__all` = ⌈10^(1.5:0.1:3.5)⌉ in 1-D and ⌈10^(1.5:0.1:2.3)⌉ in 2-D.
- θ, κ and α are `varied` keys.

Method labels are `<stepper>MUSCL<K>(<variant>,<AR|RU>,MLS<k_div>)`:
- AR is the aligned interface flux, RU the geometric one (plain Rusanov).
- MLS 0 means coupled, MLS k means frozen at k.

**The defect is a statistic, added by the post-processing pass**: `postprocess = :conservation_defect!`, loaded via
`sources = [defect_stats.jl]`. A regular `calc_stat` only sees the states of one time step, without the
parameters, so it cannot form m(t) − m(0) − B(t).

Every run file gets:
- `:conservation_defect` (time series), defined as (m(t) − m(0) − B(t)) / |m(t)|;
- `:final_conservation_defect` (scalar);
- `:boundary_flux` (time series), B(t).

Here m = h Σ u_i, the `:mass` of the hyperbolic preset. B(t) = 0 on periodic domains. For the 1-D outflow
Riemann problem, B(t) = (f(u_L) − f(u_R)) t, with f taken from the run's own equation, so Burgers, advection and
1-D Euler all work.

Check of the definition (N = 316, θ = 0.4, MUSCL2): the scheme conserves h Σ u (defect 8e-4). With the Shepard
volumes, Σ V_i u_i gives −2.9e-2, because Σ V_i = 0.955 L. So the pipeline's `:mass` is the right quantity.

"Fitted rate" below is the least-squares slope of −log|d| over log N for N ≥ 100. Values are the defect at
N = 3163 (1-D) or N = 200 per side (2-D).

## Results vs. the paper's text

| Paper statement | This run | |
|---|---|---|
| Lattice: zero to round-off or exactly −h u_L/(2m) (−2.0e-4 at 3163), the same for all schemes, decreasing like h | All 15 schemes are identical for N ≥ 80: exactly 0 at N = 100, 159, 200, 631–1259 and exactly −h u_L/(2m) elsewhere (−1.98e-4 at 3163); below N = 80 they scatter (0.8–1.5×) | ✓ |
| Frozen k_div = 2: −3.9e-4 … −6.7e-4 for K = 2–5, both fluxes, rates 0.7–0.9 | −3.9e-4 … −6.7e-4 (AR = RU to 2 digits), rates 0.71–0.90 | ✓ |
| Coupled, geometric flux: K = 3, 4 change the mass by 2.7–3.6 %, flat | −2.5 / +3.5 % (K = 3), +2.7 / +3.2 % (K = 4) (EPD1 / EPD2), flat (rate ≈ 0) | ✓ (2.5–3.5 %; the 2.7–3.6 % came from the classical-RK4 runs of 10-24) |
| Coupled, aligned flux: 6–16 %, K = 5 with EPD1 loses all its mass | 5.8–16 % for K = 3, 4, K = 5 EPD2 +15 %; **K = 5 EPD1: +20 % up to N = 1996, then the runs crash at N = 2512 and 3163 (stop at t ≈ 0.4)** | ✓ (the K = 5 EPD1 wording needs updating) |
| Frozen k_div = K, aligned flux: K = 3, and K = 4, 5 with EPD2, between 8e-5 and 9e-4 | K = 3 EPD1 −8.1e-5 ✓; **EPD2: +1.35e-3 (K = 3), +1.29e-3 (K = 4), −1.50e-3 (K = 5)** | **✗, see below** |
| Frozen k_div = K, EPD1: 1.5e-3 (K = 4), 7e-3 (K = 5), not decreasing | −1.49e-3 (rate 0.27), −6.98e-3 (rate 0.12) | ✓ |
| Unstabilized and limiters at K = 5: 2–4e-3, flat | −2.3e-3 (plain), −3.6e-3 (VK), −3.8e-3 (minmod) | ✓ |
| θ = 0.1 … 0.5: unstabilized K = 2, 4 between 1.0e-4 and 8.0e-4; MUSCL2(EPD1) −5.7e-4 … 1.7e-4 | 1.0e-4 … 8.0e-4; −5.7e-4 … +1.7e-4 | ✓ |
| MUSCL4(EPD1, k_div 4): ≤ 3.1e-4 for θ ≤ 0.3, −1.5e-3 / −1.9e-3 at θ = 0.4 / 0.5 | ≤ 3.1e-4; −1.49e-3 / −1.85e-3 | ✓ |
| κ: largest at 3.5 (3.1e-3), smallest at 9.5 (5.3e-4); geometric flux unstable at κ = 3.5 for N ≥ 1996 | +3.05e-3 / −5.3e-4; RU κ = 3.5: −1.8e-2 up to N = 1585, then u_min = −0.24 and defect −0.26 → −1.0 from N = 1996 | ✓ |
| α ∈ [0.5, 4] changes it by up to 1.3e-3 | −1.49e-3 (α = 1) … −1.8e-4 (α = 4), spread 1.3e-3 | ✓ |
| Smooth 1-D: ≤ 7.4e-8 (K = 2), ≤ 1e-9 (K ≥ 3) at 3163, orders 2–12 | 7.4e-8 / 9.6e-10 / 2.3e-12 / 4.9e-14 (θ = 0.4); fitted rates 3.3–6.8 | ✓ (the paper's 2–12 are local orders, a noisier estimator) |
| Smooth 2-D: lattice exact, ~1e-4 at 200² ∝ θ, independent of K and stabilization, orders 1.1–2.5 | lattice ≤ 3e-16; θ = 0.2: 0.8–1.4e-4, θ = 0.4: 2.1–2.6e-4 for all six methods; fitted rates 0.9–2.0 | ✓ |

### The one disagreement: frozen k_div = K with EPD2 and the aligned flux

The paper's numbers for this case are those of the **geometric** flux:

| | AR (this run) | RU (this run) | AR, old run of 10-27 | RU, old run of 10-27 | 10-27 table / paper |
|---|---|---|---|---|---|
| K = 3 | +1.35e-3 | +8.68e-4 | +1.347e-3 | +8.682e-4 | +8.7e-4 |
| K = 4 | +1.29e-3 | +6.07e-4 | +1.295e-3 | +6.067e-4 | +6.1e-4 |
| K = 5 | −1.50e-3 | −3.64e-4 | −1.497e-3 | −3.643e-4 | −3.6e-4 |

The old runs are reproduced to 4 digits (read directly from `data/review/mass/data`). The error is in the old
tabulation: `tools/mass_table.jl` groups runs by a `Key` that has no flux field. So the RU file overwrote the AR
file of the same (K, variant, MLS). This affects only the cells where the 10-26 queue had both fluxes:
EPD2, MLS = K, K = 3–5. (EPD2 MLS 2 at K = 2 agrees between the fluxes.)

Corrected statement:
- With k_div = K and the aligned flux, EPD2 leaves 1.3–1.5e-3, flat (rates 0.14–0.94), the same level as EPD1 at
  K = 4.
- Only K = 3 with EPD1 (−8e-5) reaches the floor.
- The geometric flux with EPD2 gives 3.6e-4 … 8.7e-4.

`tools/mass_table.jl` is not fixed in this branch. Its `Key` needs the flux if it is used again; the new
pipeline doesn't use it.

## Checks

- Finished: 2374/2374 runs stored with the defect statistic, 0 missing. The process exited normally, with no
  errors in the log (only a Julia world-age deprecation warning).
- Two runs did not reach t_end (MUSCL5 EPD1/AR/coupled at N = 2512 and 3163). These are the unstable
  combination, not tooling failures; they appear as "failed" in the tables.
- RU κ = 3.5 runs from N = 1996 on reach t_end but are unstable (u_min ≈ −0.24), as in the paper.
- Reproducibility: the 10-27 runs of the same configuration are matched to 4 digits; parallel and serial
  execution give identical results.

**Verdict: PASS.** The paper's numbers are reproduced through the pipeline, except the aligned-flux EPD2
frozen-k_div = K values. Those were mislabelled geometric-flux values in the old table and need correcting in
the text (1.3–1.5e-3), and the wording for coupled MUSCL5(EPD1) needs updating.
