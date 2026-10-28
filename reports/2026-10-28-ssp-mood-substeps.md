# 2026-10-28: Proof review — MOOD with SSPRK(5,4) and the first-order orientation

Review of the three implementation points of the proof check. Two are implemented in HyperCloudCore branch
**`claude/ssp-mood-substeps`** (on main `006b7fd`, which has the u2 rework); the third is wording. All runs local,
Float64, 8 threads, jitter 0.2, κ = 2.5, MUSCL K + AlignedRusanov, frozen MLS 2, EPD2/U2 (δ = 10⁻⁴Δx^D, plateau
10⁻³), scripts in `tools/ssp_mood/`.

## 1. Stage-wise MOOD and SSPRK(5,4): confirmed, implemented

The canonical Shu–Osher form (Ketcheson 2008) of each explicit tableau, `tools/ssp_mood/ssp_check.jl`:

| Tableau | SSP coefficient C | stages as convex combinations of Uⁿ and W_j = Y_j + (Δt/C) K_j |
|---|---|---|
| Euler, SSP-RK2, SSP-RK3 | 1 | every W_j enters exactly one stage |
| Ralston RK2 | 1/2 | (SSP, with half the step) |
| SSPRK(5,4) | 1.50818 | Uⁿ⁺¹ = 0.2380 Uⁿ + 0.2874 W₂ + 0.1337 W₄ + 0.3409 W₅ |
| classical RK4, IMEX pairs | 0 | — |

The final update of SSPRK(5,4) reuses the substeps of stages 2 and 4, which were only checked through
Y₃ = 0.444 Uⁿ + 0.556 W₂ and Y₅ = 0.178 Uⁿ + 0.822 W₄; an admissible Y₃ does not make W₂ admissible, so the
convex-combination argument of Corollary 6.1 does not cover it. (The reviewer's decomposition has a negative Uⁿ
coefficient; the canonical one has none, but the gap is the same: the substeps entering the final update are not
checked.) For SSP-RK2/RK3 the stage check is sufficient for admissibility, but the substep is also the natural
object of a DMP check against Y_{k−1} (Y₃ = ¾Uⁿ + ¼W₂ of RK3 mixes in the old state).

**Change**: for an explicit tableau with C > 0 MOOD now checks W_{k−1} = Y_{k−1} + (Δt/C) K_{E,k−1} at checkpoint k
(against the neighbourhood of Y_{k−1}, PAD and criterion), instead of Y_k; the recompute set is unchanged (only
K_{k−1}). Then every stage and the final update are convex combinations of Uⁿ and checked states, for any SSP
method. `ssp_coefficient(bt)` computes C by bisection on the canonical form; non-SSP and IMEX tableaus keep the
candidate check; `:Time_mood_substeps = false` restores it everywhere. Cost: one axpy kernel per check, one more
state per particle when MOOD is on.

Results (off → on):

| Case | K = 2 (RK2) | K = 3 (RK3) | K = 4 (SSPRK54) |
|---|---|---|---|
| Sod 400, L1(ρ) | 6.489e-3 → 6.487e-3 | 6.308e-3 → 6.409e-3 | 5.311e-3 → 6.059e-3 |
| 1-2-3 problem 400, L1(ρ) / min p | 8.16e-3 / 6.04e-3 → 7.65e-3 / 6.46e-3 | 5.42e-3 / 6.47e-3 → 5.35e-3 / 6.60e-3 | 5.87e-3 / 6.81e-3 → 5.81e-3 / 6.82e-3 |
| near vacuum (u = ∓3.5), CFL 0.4, min p | – | 4.51e-3 → 4.52e-3 | 4.58e-3 → 4.59e-3 |
| near vacuum, CFL 0.8, min p (L1) | – | 4.33e-3 → 3.36e-3 (1.03e-2 → 1.12e-2) | 4.54e-3 → 4.63e-3 (7.98e-3 → 8.29e-3) |
| box 2-D 100², over / undershoot | 9.5e-5 / 2.2e-4 → 3.2e-4 / 3.1e-4 | 3.5e-4 / 4.5e-4 → 3.8e-4 / 4.6e-4 | 3.0e-4 / 3.4e-4 → 3.3e-4 / 3.2e-4 |
| vortex 141² (plateau 10⁻²), L1(ρ) | – | 0.1049 → 0.1049 | 0.0603 → 0.0603 |

Runtimes are unchanged. No case produced an inadmissible state with either check (every run positive, also near
vacuum): in these tests the gap is a gap of the guarantee, not an observed failure. The substep check is the more
conservative one (it checks the full Euler step), which costs up to 14 % in the Sod error of SSPRK(5,4) and changes
the others by about ±1–10 %. Kept on by default because it makes the PAD/DMP statement hold for every SSP tableau.

**For the paper**: state the result for the substeps (the reviewer's suggested proposition: particle i at order 1,
EPD1/EPD2/StrictEPD0, frozen mode with the aligned flux, Y_{k−1} admissible on ν̄(i) and (5.5) at Y_{k−1} ⇒ the
substep at i is admissible), and note that Corollary 6.1 as stated covers only the all-order-one case.

## 2. Aligned orientation at K = 1: confirmed, implemented

MUSCL applied the aligned orientation only on edges where a particle is below the full order, which never happens
for K = 1. On a jittered 2-D cloud (24², jitter 0.2, κ = 2.5, degree-1 weights) the geometric orientation gives 200 of
11 258 first-order coefficients a_ij < 0 (1.8 %, down to −0.57); with the aligned orientation none. **Change**: the
aligned orientation is also used on every edge of a first-order scheme. Unit test (fails on main, passes now).
Effect on runs (RK2MUSCL1): Sod and the 1-2-3 problem are bit-identical (1-D); the 2-D box L1 1.0729e-1 → 1.0739e-1,
bounds unchanged (both runs bounded). Like for the substeps, the change closes a gap of the guarantee (Proposition 5.2
for K = 1) that the runs did not exhibit.

Corollary 5.3: the "Hence …" sentence should be restricted to particles with k_i < K (under EPD2 an order-K particle
next to a lowered one keeps the geometric orientation on its edges to other order-K particles).

## 3. Where the aligned signs come from: wording

The signs are read from the cached derivative weights of the divergence fit in use: with a fixed `div_order` (frozen
mode) that is the degree q(k_div) fit, i.e. degree 1 for k_div = 2, as Section 7.4 says. Only with `div_order = 0`
(coupled mode) does the cache belong to the full degree q(K) while lowered particles fit with a lower degree, and then
condition (b) is not guaranteed. Remark 5.4(v) should say "the derivative weights of the divergence fit in use (degree
q(k_div) in frozen mode)". HyperCloudCore's docs now say so.

## Tests

HyperCloudCore: full suite on the branch **787/787**. New tests:
- SSP coefficients and the canonical form;
- the recorded MOOD candidates are the substeps (RK3, SSPRK54), or the stage candidates with `mood_substeps = false`;
- the K = 1 coefficients (this test fails on main).

HyperCloud suite against the branch (local dev env, no GUI): 920/924. The 4 failures are the `create_method` label
tests, which fail on HyperCloud main since b061464 (default flux `:AlignedRusanov`), unrelated.

## Verdict

**PASS**: both implementation points are correct and implemented; neither degrades the runs beyond the stricter
substep check of SSPRK(5,4) on Sod (+14 % L1); no positivity failure was observed with the old checks either.
