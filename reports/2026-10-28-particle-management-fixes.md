# 2026-10-28: Particle management, review of the proposals and the HyperCloudCore fixes

Follow-up to `reports/2026-10-28-particle-management.md` (PR #11). The five proposed HyperCloudCore changes were
reviewed against the code; four are implemented in HyperCloudCore branch **`claude/adaptivity-conservation`**
(commit `fad1503`, on top of main `2fe8b46`). The study was rerun locally with the same driver
(`tools/adapt_study.jl`, unchanged) against main and the fix branch: 47 runs, all completed, Float64,
8 threads, `tools/adapt_fix_batch.sh`. Raw RESULT lines: `results/adapt_fix/all.log` (gitignored); tables:
`reports/2026-10-28-adapt-fix-tables.md` (`python3 tools/adapt_fix_table.py results/adapt_fix/all.log`).
The main-branch numbers reproduce the cluster study (e.g. lin, adaptivity on, N = 1600: −4.18e-5, report −4.2e-5;
b1pre N = 200: re-volume part −6.34e-4, report −6.3e-4).

## Review of the proposals

| # | Proposal | Verdict | Implementation |
|---|---|---|---|
| 1 | Pass the grid's volume estimator to the stepper | **accepted**, confirmed: with `Grid_volumes = :mls` the stepper of the run pipeline re-estimated with `ShepardVolumes` | `build_particle_grid` stores `context[:VolumeEstimator]`; stepper and adaptivity builder read it |
| 2 | Event spacing over all changed particles (periodic) | **accepted**, cause found (below) | spacing `R (2 + 2√D) + dx` instead of `4R`; `_collect_near!` already wraps periodic axes |
| 3 | Consistent volumes between remap and estimator | **accepted**, in a different form than "keep the content" (below) | `ParticleAdaptivity(volume_estimator = …)`, passed by the builder on moving grids |
| 4 | Expose `max_events` | **accepted** | `:Adapt_max_events` (`build_adaptivity`) |
| 5 | MLS volumes as default on moving grids | **not accepted as a default** (cost, device path) | documented the Shepard bias and the MLS cost; MLS on a device grid now works through a host copy |

**Cause of #2.** The events of one call are remapped one after the other on the *old* topology; a merged-away particle
is only marked and stays in the cloud until compaction. The local set of a remap (every particle whose shape
function can change) reaches `R (1 + 2√D) + dx` from its centre, i.e. 3.4R in 1-D but 4.2R in 2-D, while
accepted centres could be as close as `4R − R/2 = 3.5R`. In 2-D a later remap therefore gave part of its swept
content to a particle an earlier event had removed, and that share disappeared at compaction. The new spacing
keeps every event's remap region free of the particles another event removes or adds. A unit test with close
pairs 4.2R apart (2-D): worst call **1.75e-7 → 1.3e-15** of the mass (main was only exact with `max_events = 1`).

**#3, why not keep each particle's content.** The remap's volumes are integrals of the Shepard partition; the
estimators define volumes differently (MLS degree 2, or the particle-based Shepard KDE). Next to a merge the two
differ by up to 18 % (MLS) / 9 % (Shepard) per particle, so `U_i ← V_i^remap U_i / V_i^est` is conservative but
puts that difference into the states (consistent with the NaNs `ADAPT_KEEP_CONTENT` produced). Implemented
instead: right after the topology change the volumes are re-estimated with the stepper's estimator, and the
resulting change of `Σ V U` is put back onto the particles the remap changed with the remap's own bounded,
exactly conservative repair. In the unit test the states then move by 3.2e-4 for a defect of 2.7e-4 (instead of 18 %).

**#5, cost.** Per-step MLS volumes (degree 2) on a Lagrangian Burgers grid, 8 threads: 3.0× a MUSCL2/RK2 step in
1-D (N = 1600) and **21×** in 2-D (141²), against 5 % for Shepard; the MLS sweep is a serial host loop and was not
usable on device grids. Making it the default would make 2-D moving-grid runs ~20× slower. A kernelised MLS
estimate (two passes over the background lattice) would be the way to make it default-worthy.

## Results: main vs fix (MLS volumes, K = 2, relative mass change at t_end)

"Re-volume part" is what re-estimating the volumes after the adaptation calls changes (the driver's
`rel_revolume_sum`); "max call" the largest mass change of one adaptation call (absolute).

| Case | N | adapt off | on, main | on, fix | re-volume part main → fix | max call main → fix |
|---|--:|--:|--:|--:|--:|--:|
| b1pre (1-D Burgers, smooth) | 400 | −1.8e-6 | +5.9e-5 | −5.5e-6 | +8.3e-5 → 0 | 9e-16 → 1e-15 |
| | 800 | +9.9e-8 | +1.4e-5 | −2.1e-6 | +2.2e-5 → 0 | 2e-15 → 1e-15 |
| | 1600 | +4.8e-8 | −4.4e-6 | **−1.6e-7** | −4.4e-6 → 0 | 2e-15 → 2e-15 |
| lin (1-D, prescribed motion) | 400 | −4.0e-4 | −6.5e-4 | −4.1e-4 | −2.3e-4 → 0 | 1e-15 → 2e-15 |
| | 800 | −6.1e-5 | −1.3e-4 | −1.1e-4 | −1.8e-5 → 0 | 2e-15 → 2e-15 |
| | 1600 | −4.8e-6 | −4.2e-5 | −2.6e-5 | −1.6e-5 → 0 | 2e-15 → 3e-15 |
| vortex (2-D Euler) | 71 | +8.0e-4 | +7.8e-4 | +8.5e-4 | −2.7e-5 → 0 | **2.5e-5 → 1.4e-12** |
| b1 (1-D Burgers, after shock) | 200 | −9.4e-3 | −6.0e-3 | −1.3e-2 | +5.0e-3 → 0 | 1e-15 → 2e-15 |
| | 800 | −1.6e-2 | −1.1e-2 | −1.2e-2 | −5.3e-4 → 0 | 2e-15 → 2e-15 |
| b2 (2-D Burgers, after shock) | 100 | +4.5e-3 | +1.7e-3 | +1.1e-2 | −8.1e-3 → 0 | **2.3e-10 → 5.3e-14** |

- **Every adaptation call is exact** with the fix: max call ≤ 5.3e-14 absolute in the Burgers and advection runs and
  ≤ 1.5e-12 absolute in the vortex runs (total mass about 100, i.e. 1.5e-14 relative), and the **re-volume part is
  exactly 0** in every run.
- **Smooth flows**: on b1pre adaptivity no longer costs 1–2 orders: at N = 1600 −4.4e-6 → −1.6e-7 (off: +4.8e-8),
  converging again (fix: −5.5e-6, −2.1e-6, −1.6e-7). On lin the volume part is gone (−4.2e-5 → −2.6e-5 at 1600),
  but adaptivity still costs about a factor 5 over the unadapted run (−4.8e-6). That remainder is no longer the
  remap (exact) but the scheme on the cloud with inserted particles: `scheme_vs_exact_rel` is −2.6e-5 with vs −4.8e-6
  without adaptivity, and it converges at second order (4.1e-4, 1.1e-4, 2.6e-5), the remap's accuracy order,
  while the unadapted run converges faster.
- **After shocks** the fix is not better and sometimes looks worse (b1 N = 200: −6.0e-3 → −1.3e-2; b2 100²:
  +1.7e-3 → +1.1e-2). On main the re-volume error had the opposite sign of the Lagrangian shock loss and partly
  cancelled it (b1 N = 200: +5.0e-3, b2: −8.1e-3). With it removed, adaptivity is roughly neutral after shocks in
  1-D (on vs off: −1.3e-2 / −9.4e-3, −7.9e-3 / −1.1e-2, −1.2e-2 / −1.6e-2) and does not help in 2-D. **The earlier
  finding "adaptivity halves the shock loss in 1-D" was partly this cancellation.** The percent-level, non-converging
  loss after shocks is the ALE update itself (adaptivity off shows it too), as the original report said; nothing in
  this change addresses it.
- **Cost**: Shepard runs +2 % (vortex 71²); MLS runs +5 % (1-D) to +22 % (vortex 71²), the extra MLS estimate per
  adaptation call. The wider spacing allows fewer events per call in 2-D, so the cloud refines more slowly
  (b2 100² by t = 1.5: 894 instead of 1402 insertions). All runs finish. On 2-D Burgers with adaptivity main
  went slightly negative (min u₁ −2.2e-5 at 71², −4.0e-7 at 100²); with the fix all Burgers runs stay positive
  (b2: +7.0e-6, +1.7e-4; 1-D: 1.4e-11 as without adaptivity).

## Not changed

- Defaults: Shepard volumes stay the default (their resolution-independent bias on moving grids is now documented
  in HyperCloudCore; use `Grid_volumes = :mls` when the mass of a moving grid matters).
- The shock band exclusion (`exclude_reduced_order`) still uses 4R; the remap's reach in 2-D is larger, but this
  only affects which events are allowed, not conservation.
- The driver (`tools/adapt_study.jl`) still works around bug #1 and builds `ParticleAdaptivity` itself for
  `ADAPT_MAXEVENTS` (that path drops the volume estimator with the fix branch; use `Adapt_max_events` instead).

## Verdict

The four accepted library fixes do what they should: **PASS** for exact conservation of every adaptation call and
across the following volume update, and for the restored convergence of the mass on smooth adapted flows.
**INCONCLUSIVE** for shocked Lagrangian flows: the mass loss there is the moving-grid update, unchanged.
