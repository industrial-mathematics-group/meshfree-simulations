# 2026-10-28: Why MOOD/u2 drops the isentropic vortex to first order, and the fix

Follow-up to the review rerun (`reports/2026-10-27-review-rerun.md`, `2026-10-27-vortex-review-table.md`): with
MOOD (EPD2, U2) every scheme was of first order on the vortex, errors 100–4·10⁴× those without MOOD. Fix in
HyperCloudCore branch **`claude/mood-u2-extrema`** (commit `e92b940`, on main `463172b`). All runs local, Float64,
8 threads, the review setup (β = 5, [0, 10]², jitter 0.2, κ = 2.5, CFL 0.8, t = 10, MUSCL K + AlignedRusanov,
frozen MLS 2, EPD2/U2 with `MOOD_delta_relax = 1e-4`, i.e. δ = 10⁻⁴Δx²). Scripts: `tools/mood_u2/` (README there).
The main-branch numbers reproduce the review table (K = 3, N = 100: 2.015, review 2.024 in Float32).

## 1. Diagnosis

An instrumented copy of the u2 criterion (`tools/mood_u2/u2_override.jl`) counts every DMP violation and why u2
rejected it. N = 100, K = 3, first steps (t = 0.1):

- **The violations are not at the core.** 13 903 of 14 264 are in the density, at distance **2.75–4.5** from the
  centre (82 within 0.5). There ρ − 1 ≈ 0.23 e^{1−r²} is 8·10⁻⁵ … 7·10⁻⁸: the density varies by 1–1000 δ across a
  stencil, and the candidates overshoot by 1–10 δ (δ = 10⁻⁶ at N = 100).
- **Why the density overshoots there.** Its flux is the momentum, and the velocity perturbation decays only like
  e^{−r²/2} (4·10⁻² at r = 3): the density update carries the truncation error of div(ρu), ≈ 10⁻⁶, in a region where
  ρ itself is nearly flat. The DMP on ρ cannot tell that noise from a new extremum.
- **u2 rejects them all** (6 of 14 264 accepted): the sign test rejects (curvature ∂ₓₓρ of the tail passes through
  zero along the lines |x − x_c| ≈ 0.7), the ratio test |c|min/|c|max ≥ ½ rejects (over a stencil of radius 2.5h the
  tail's curvature decays by e^{−3} at N = 100), curvature was tested also in components that did not violate the
  DMP, and once a particle is at order 2 u2 is unavailable, so it falls to order 1 (8 959 rejections). Within t = 0.1
  the density error is already 0.28 and the core minimum has risen from 0.49 to 0.52.

The paper reviewer's two remarks are both confirmed: (i) the flat-region guard δ³ of (5.3) is ≈ 10⁻¹⁸ (N = 100) to
2·10⁻²² (N = 400) and never applies, so the relaxed DMP is essentially strict; (ii) (5.4) compares δ with curvatures
and with products of curvatures (mixed units). Fixing (ii) alone changes nothing measurable (2.015 → 2.015).

## 2. What did and did not work

Prototype variants (N = 100, K = 3; vortex density L1 at t = 10, no MOOD 0.229):

| Variant | Vortex | Box (2-D, K = 3, 100²) over/undershoot |
|---|---|---|
| main | 2.02 | 5e-5 / 4e-7 |
| units fixed | 2.02 | |
| + only violating components | 2.01 | |
| + no ratio test | 1.33 | |
| absolute δ = 10⁻⁴ (+ the above) | **0.229** (= no MOOD) | 5.6e-3 / 1.1e-2 (each stage may exceed the current bounds by δ: the overshoot ratchets) |
| plateau test η = 10⁻³ / 10⁻² of the global range, no ratio test | 1.30 / **0.229** | 4.0e-4 / 6.8e-4 ; 1.3e-2 / 8.1e-3 |

A larger absolute δ is not the answer (discontinuities ratchet). A **plateau test** relative to the component's
global range keeps discontinuities bounded (where a component jumps its local range is the jump) and accepts the
tail. Dropping the ratio test entirely is wrong: on a 1-D lattice a box then rings by **9.7 %** (HyperCloudCore's
pipeline boundedness test caught it): the plateau next to a jump has exactly zero curvature, so the sign test passes
and only the ratio test sees the collapse. It is kept, but only among the nearest neighbours (1.5h), the face
neighbours it was designed for.

## 3. The fix (`MOODu2(δ; plateau = 1e-3, ratio = 1/2)`, keys `:MOOD_plateau`, `:MOOD_curvature_ratio`)

1. **Plateau**: a component that varies across the stencil by less than `plateau` × its global range is accepted
   (one reduction per MOOD check), also below order 3.
2. Curvature tests only on the components that **violate the DMP**.
3. Curvatures compared as the **variation c·h² over one spacing** (sign: (c_min h²)(c_max h²) > −δ²; negligible:
   |c|max h² < δ).
4. Diot's **ratio test among the nearest neighbours** (≤ 1.5h); sign test over the whole stencil.

## 4. Results with the defaults

Vortex, density L1 vs IC at t = 10 (min ρ in brackets; exact 0.49):

| N | K = 3 no MOOD | K = 3 main | K = 3 fix | K = 3 plateau 1e-2 | K = 4 no MOOD | K = 4 main | K = 4 fix | K = 4 plateau 1e-2 |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 50 | 0.949 | 2.658 | 2.625 | 2.611 | 0.419 | 2.539 | 2.488 | 2.465 |
| 71 | 0.483 | – | 2.309 | 2.217 | 0.229 | – | 2.176 | 2.100 |
| 100 | 0.229 | 2.015 | 1.912 | 1.033 | 0.119 | 1.863 | 1.761 | 0.857 |
| 141 | 0.105 | – | 0.942 | **0.105** | 0.060 | – | 0.953 | **0.060** |
| 200 | 0.0455 (0.490) | 1.218 (0.730) | **0.0455 (0.490)** | – | 0.0303 (0.494) | 1.119 (0.708) | **0.0303 (0.494)** | – |

From 200² (default) or 141² (plateau 10⁻²) MOOD leaves the vortex untouched: the errors equal the no-MOOD errors to
all printed digits. On the coarser clouds (≤ 100², ≤ 10 particles across the core radius) the nearest-neighbour
curvature test still rejects particles, and MOOD still costs accuracy there.

Discontinuities (library, defaults vs main):

| Case | main over/undershoot | fix over/undershoot | L1 main → fix |
|---|---|---|---|
| 1-D box, lattice, N = 100, K = 3 / 4 / 5 (Core's test, tolerance 1e-2); left column: fix with `plateau = 0` | 5.9e-5 / 1.6e-5 / 2.5e-6 | 1.6e-3 / 6.0e-5 / 1.8e-4 | |
| 2-D box, jitter 0.2, 200², K = 3 | 1.6e-7 / 1.5e-7 | 1.1e-3 / 8.7e-4 | 1.88e-2 → 1.88e-2 |
| 2-D box, 200², K = 4 | 1.9e-7 / 1.5e-7 | 6.2e-4 / 5.7e-4 | 1.47e-2 → 1.48e-2 |
| Sod, 400, jitter 0.2, K = 3 | 4e-7 / 1e-7 | 2.5e-4 / 9e-6 | 6.29e-3 → 6.31e-3 |
| Sod, 400, K = 4 | 7e-14 / 8e-9 | 1e-13 / 2.7e-6 | 5.89e-3 → 5.31e-3 |

Discontinuities stay bounded to about 10⁻³ of the jump (before 10⁻⁷); L1 errors unchanged or better.

Library tests: HyperCloudCore 732/732 (16 new u2 tests); HyperCloud 920/924 against the branch (the 4 failures are
`create_method` label tests that expect the old default flux; they fail on HyperCloud main since commit b061464,
unrelated to MOOD).

## 5. For the paper

- Section 5.3/5.4: the DMP needs a flat-region criterion that actually triggers (the plateau test, relative to the
  global range of the component), and u2's curvature tests need consistent units (c·h² against δ). The decisive
  changes for the vortex are the plateau test, the restriction to the violating components and the stencil-adapted
  ratio test; the unit fix alone does not change the result.
- Section 7.3: with the fix, MOOD keeps the vortex at full order from N = 200 (141 with plateau 10⁻²); the
  remaining coarse-grid degradation is the curvature test on under-resolved extrema. The background velocity in
  the experiment file is (1, 1) (the TODO in the text).
- Results computed with u2 before this change are not reproduced by the new default (hashes of the parameter dicts
  are unchanged); rerun what goes into the paper.

## Verdict

**PASS** for the fix at the resolutions where the vortex is resolved (N ≥ 200 default, N ≥ 141 with plateau 10⁻²:
errors identical to no MOOD) with discontinuities bounded to ~10⁻³; **INCONCLUSIVE** for N ≤ 100, where MOOD still
degrades the vortex. Not run: Float32, N = 283/400 (too long locally), GPU.
