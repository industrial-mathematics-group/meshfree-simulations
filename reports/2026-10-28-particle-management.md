# 2026-10-28: Particle management (insertion/merging) and mass conservation

HyperCloud 0.5.0 / HyperCloudCore 0.7.3, CPU, Float64. 70 cluster jobs (778 runs, all completed; queue
`cluster/queues/2026-10-28-adapt.txt`; the RESULT lines of all jobs in `reports/2026-10-28-adapt-results.txt`) plus local diagnostic runs.
Full per-N tables: `reports/2026-10-28-adapt-tables.md` (`python3 tools/adapt_table.py <logs>`).

**Tooling.** `tools/adapt_study.jl` builds a run exactly like `run_direct_simulation` but with its own time loop, so every
`adapt_particles!` call is recorded with its own `mass_residual` (the change of Σ V_i U_i over the call), separately
from the change over the time steps. Diagnostic switches: `ADAPT_EVENTLOG` (events of non-conservative calls),
`ADAPT_HOLECHECK`, `ADAPT_MAXEVENTS`, `ADAPT_KEEP_CONTENT`. For linear advection it also evaluates the discrete mass of the
*exact* solution on the final cloud (`quad_exact_rel`), which separates the scheme from the measurement.
`cluster/run_adapt.sh` runs one case / order / volume estimator / adaptivity setting over a range of N.

Cases (all periodic, so the exact mass is constant; jitter 0.2; R = 2.5): **lin** 1-D advection of a Gaussian with a
prescribed sinusoidal grid velocity (A = 0.5, t = 2); **b1pre / b1** 1-D Burgers Gaussian with the Lagrangian mover,
t = 0.5 (smooth) / t = 2 (after the shock); **b2** 2-D Burgers Gaussian, Lagrangian, t = 1.5 (after the shock);
**vortex** isentropic vortex, Lagrangian, one period. Schemes: RK2MUSCL2 / RK4SSPMUSCL4 with EPD2, AlignedRusanov,
frozen MLS 2 (the recommended default). Adaptivity defaults: max_gap 1.6, min_gap 0.5, samples_per_dx 2,
exclude_reduced_order on.

## 1. The remap itself: exact per event, but simultaneous events interact

- Every insertion/merge call in 1-D, and almost every call in 2-D, changes Σ V_i U_i by ≤ 1e-14 (relative), for every
  component of the Euler state as well: the swept-region remap is conservative by construction, as documented.
- **But individual calls with several events at once are not**: up to 2.5e-5 absolute (vortex, N = 71; accumulated
  7e-8 … 8e-7 relative) and 1e-8 (2-D Burgers, N = 283). The bad call on the vortex had four events, two of them
  4.1 R apart across the periodic seam; the cloud around them was fully covered (no point of the quadrature box
  without a particle within R), so it is not a partition-of-unity hole. **With `max_events = 1` the same run is
  exact again** (accumulated remap residual 1e-14, max per call 6e-13). So the "events are kept 4R apart" rule does
  not fully decouple the remaps; my guess (not verified in the code): the rule measures from the event *center*,
  while the removed partner of a merge (or the shape functions of the particles between two events) reaches into
  the other event's region. Fix in HyperCloudCore: apply the separation to all changed particles (periodic-aware)
  with a margin, or process events whose influence regions overlap in separate passes. Stop-gap: `max_events = 1`
  (not exposed by `build_adaptivity`).

## 2. The remap's conservation is undone by the next volume update

The remap conserves Σ V_i U_i **with the volumes it computes for the changed particles** (δV from the shape
functions). At the end of the next step the volume estimator recomputes all volumes from the geometry, and that
changes Σ V_i U_i. Measured right after each adaptation call (1-D Burgers, smooth, MLS volumes): the re-estimation
alone changes the mass by −6.3e-4 (N = 200, total defect −1.1e-3) and +2.2e-5 (N = 800, total +1.4e-5), i.e. it is
the dominant part of what adaptivity costs. Consequence (cases with MLS volumes, mass change at N = 1600):

| Case | adaptivity off | adaptivity on |
|---|---|---|
| lin (smooth, prescribed motion), K = 2 / 4 | −4.8e-6 / −1.9e-5 (order 3.9 / 2.1) | −4.2e-5 / −4.6e-5 (order 2.1 / 1.3) |
| b1pre (Burgers, smooth), K = 2 / 4 | +4.8e-8 / +5.3e-8 | −4.4e-6 / −6.5e-6 |

So on smooth flows adaptivity costs about one to two orders of magnitude in mass conservation (still converging).
Keeping each particle's *content* when the volumes are re-estimated (`ADAPT_KEEP_CONTENT=1`: U_i ← U_i V_i^remap /
V_i^new) is **not** a fix: it stops converging (1.7e-3 … 2.8e-4) and produces NaNs after the shock. A consistent fix
needs the remap and the volume estimator to agree, e.g. the remap producing the estimator's volumes for the changed
particles (and their neighbours), or the estimator not overwriting volumes of untouched particles.

## 3. Moving grids: Shepard volumes are not a convergent quadrature (and a wiring bug)

- On a lattice moved by a smooth prescribed velocity (no MOOD, no adaptivity), Σ V_i u_i of the **exact** solution
  with Shepard volumes is off by **−2.0e-3 for every N (200 … 1600) and every order**, also with 0 or 3 smoothing
  sweeps, while the scheme converges to the exact solution at second order (difference 1.0e-3 → 6.1e-5 → 1.5e-5).
  With **MLS volumes** the quadrature converges (degree 2: −2e-5 → 3e-8 → 5e-9) and the mass change becomes the
  scheme's own (8.6e-4 → 5.4e-5 → 1.3e-5). The docs' "ALE mass errors converge with the scheme's order" holds only
  with MLS volumes; with the default Shepard volumes a moving grid shows an O(1)-in-N mass drift.
- **Bug (HyperCloudCore 0.7.3):** `Grid_volumes = :mls` only sets the *initial* volumes. `build_timestepper` takes the
  estimator from `context[:VolumeEstimator]`, which nothing sets, so the stepper always re-estimates with the default
  `ShepardVolumes()` after each step. `tools/adapt_study.jl` sets it explicitly; the run pipeline (`run_direct_simulation`)
  is affected. One-line fix: set `context[:VolumeEstimator] = vol` in `build_particle_grid` (or pass the grid's
  estimator in `build_timestepper`).

Effect in the study (adaptivity off, mass change at N_max):

| Case | Shepard | MLS |
|---|---|---|
| lin, K = 2 | −5.1e-4 (order 0.3) | −4.8e-6 (order 3.9) |
| b1pre, K = 2 | +2.3e-3 (0.4) | +4.8e-8 (2.6) |
| b2, K = 2 | −3.8e-2 | −6.6e-3 |
| vortex, K = 2 | +3.7e-4 | +5.3e-4 |

## 4. After shocks, the Lagrangian cloud is not conservative (independent of adaptivity)

| Case (MLS volumes) | adaptivity off | adaptivity on (defaults) | best setting |
|---|---|---|---|
| b1 (1-D Burgers after the shock), K = 2 | −2.2e-2 (no convergence) | −1.0e-2 | −5.1e-3 (exclude_reduced_order off), −5.9e-3 (min_gap 0.7) |
| b1, K = 4 | −1.9e-2 | −6.3e-3 | |
| b2 (2-D Burgers after the shock), K = 2 | −6.6e-3 | −6.4e-3 | −1.5e-3 (samples_per_dx 4, at 141²) |
| b2, K = 4 | −1.1e-2 | −1.4e-2 | |
| vortex, K = 2 / 4 | +5.3e-4 / +5.9e-4 | +4.9e-4 / +4.2e-4 | |

- With a Lagrangian mover the particles pile up at the shock, and the mass error is percent-level and does not
  converge with N — with or without adaptivity, i.e. it is the ALE update, not the remap. The ALE form uses each
  particle's own grid velocity (−∇·(F − w_i U) at x_i); across a shock w jumps, so the w_i U term is no longer a
  pairwise-antisymmetric flux (my interpretation, consistent with smooth cases converging and shocked ones not).
- Adaptivity helps in 1-D (halves the error: insertions keep the rarefaction resolved, merges thin the pile-up),
  more so when it is allowed in the shock band (`exclude_reduced_order` off: −5.1e-3) or merges earlier
  (min_gap 0.7). In 2-D it is neutral for K = 2 and slightly harmful for K = 4.
- All runs stay positive (min u₁ ≥ −5e-7 on Burgers) and finish; the vortex keeps about 10–15 % more particles with
  adaptivity (40 000 → 46 770 at 200²; 23 500 insertions, 16 800 merges).

## 5. Summary

| Question | Answer |
|---|---|
| Is a single insertion/merge conservative? | **Yes**, to round-off (1e-14), in 1-D and 2-D, for all Euler components |
| Is an adaptation call with several events conservative? | **No**: interacting events, up to 1e-5 per call (2-D); exact with `max_events = 1` → HyperCloudCore fix |
| Does the conservation survive the next step? | **No**: the volume re-estimation changes Σ V U by about as much as the whole adaptivity defect → library fix needed |
| Can mass be measured on moving grids? | Only with MLS volumes; Shepard volumes have an O(1)-in-N bias (−2e-3 here) |
| Does `Grid_volumes` work? | Only for the initial volumes (wiring bug) |
| Does adaptivity help? | Smooth flows: costs 1–2 orders in mass error. After shocks: halves the (non-converging, percent-level) error of the Lagrangian cloud in 1-D, neutral in 2-D |
| Recommended now | MLS volumes (with the wiring fix), `max_events = 1` until the spacing rule is fixed; for shocked flows the Lagrangian mover itself is the limiting factor |

Proposed HyperCloudCore changes (not made, other repo): (1) pass the grid's volume estimator to the stepper;
(2) separation rule over all changed particles (or one event per neighbourhood pass); (3) consistent volumes between
the remap and the estimator; (4) expose `max_events` in `build_adaptivity`; (5) consider MLS volumes as the default
on moving grids.
