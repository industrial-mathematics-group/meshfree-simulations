# Mass conservation of the meshfree MUSCL scheme on irregular point clouds

Study of 2026-10-17/18. Everything runs on the CPU in Float64 from `study/conservation/`
(`julia --project=env_study study/conservation/<script>.jl`); raw results are in `out/`.

- `Harness.jl` builds periodic Burgers / advection runs on [0,1]^D with the library builders, adds
  weight functions (`:shifted_exp`, `:wendland`, `:gauss`) and extracts the MLS coefficients c_ij.
- `Solver.jl` is a standalone reimplementation of the library scheme (scalar laws, MUSCL of any order,
  Rusanov, SSP-RK2/RK4) that reproduces the library right-hand side to 1e-15 for MUSCL2–4 in 1-D and
  2-D (`validate.jl`, `validate_ho.jl`), plus alternative conservative formulations.

Test problem: periodic Burgers, u0 = 1 + ½ sin(2πx) [· sin(2πy)], jittered lattice (randomness
`rf`·Δx), shock from t = 1/π. "Drift" = relative change of Σ V_i U_i (Shepard volumes).

## 0. A bug first: frozen particles at the periodic seam (fixed, HyperCloudCore PR #9)

Jittered lattice points that left a periodic axis at its lower end were classified *before* the
periodic wrap and became boundary particles, which are never updated. In 2-D, 4 % of the particles
(the left and bottom seam) were frozen. This was the largest conservation error on periodic jittered
grids: 2-D Burgers drift 4e-3…7e-3 → 5e-6…2e-4 after the fix (and the extrema now match the regular
grid). All results below are with the fix.

## 1. Where the remaining error comes from

The scheme is `dU_i/dt = −2 Σ_j c_ij · (F*_ij − F_i)` with c_ij the MLS first-derivative weights of
particle i (fit of the differences, no constant term). For any weights V:

```
d/dt Σ_i V_i U_i = −Σ_pairs (V_i c_ij + V_j c_ji) · F*_ij  +  2 Σ_i V_i (Σ_j c_ij) · F_i
                     └── pair antisymmetry defect ──┘          └── row-sum defect ──┘
```

Both vanish on a lattice (symmetry) and are O(jitter) otherwise (`defects.jl`):

| 1-D, R = 2.5, degree 1 | rf = 0.01 | 0.05 | 0.2 |
|---|---|---|---|
| antisymmetry (Shepard V) | 2.6e-3 | 1.3e-2 | 5.2e-2 |
| row sums Σ_j c_ij | 1.5e-3 | 7.6e-3 | 3.0e-2 |

Larger radii reduce both (~1/R), degree-3 fits have ~3× larger antisymmetry defects. In 1-D there are
positive weights V* that make the pair term vanish exactly; the drift measured with them is unchanged
(1.2e-3 vs 1.0e-3), so **the row-sum term dominates in 1-D**.

How it shows up (library scheme, `e1_library.jl`, `e6_misc.jl`):

- Smooth solutions: drift shrinks fast with N (MUSCL4 Burgers before the shock: 8e-4 → 4e-6 for
  N = 50 → 800). Non-conservation is a shock problem.
- After the shock it is **systematic** (monotone in time, accumulating while the shock moves) and
  **does not vanish under refinement** with limited MUSCL2: 1-D, R = 2.5, rf = 0.2, VK: 0.8e-3, 1.0e-3,
  1.0e-3, 1.0e-3, 1.2e-3 for N = 50…800; it grows like rf² (rf = 0.05/0.1/0.2/0.3: 8e-5/3e-4/1.2e-3/3.4e-3).
  First order is 10× smaller. 2-D: 1.7e-4 (R = 2.5), 8e-6 (R = 5.5).
- Accuracy: the post-shock L1 error vs a 40000-cell Godunov reference still converges at first order
  (1.1e-3 at N = 800, twice the regular-grid error).

## 2. Knobs without changing the scheme

**MOOD with a coupled MLS order is by far the largest source** (`e6_misc.jl` b; MUSCL4 MOOD EPD1, R = 5.5):

| | regular grid | rf = 0.2 |
|---|---|---|
| 1-D, MLS_order = 4 (fixed) | 6e-16 | 2.2e-4 |
| 1-D, MLS_order = 0 (coupled) | **2.3e-2** | **2.3e-2** |
| 2-D, fixed | 8e-16 | 3.2e-5 |
| 2-D, coupled | **1.8e-3** | **2.1e-3** |

So keeping a fixed global MLS order is right for conservation (by ~100×).

**Weight function and support radius** (MUSCL2 VK, rf = 0.2, drift after the shock; N = 400 in 1-D, 64² in 2-D):

| | 1-D R = 2.5 | 1-D R = 5.5 | 2-D R = 2.5 | 2-D R = 5.5 |
|---|---|---|---|---|
| exponential α = 1 (library default) | 1.0e-3 | 3.4e-4 | 1.7e-4 | 1.7e-6 |
| exponential α = 2 | **5.1e-5** | 8.4e-5 | 7.2e-5 | 7.0e-6 |
| exponential α = 4 | **5.1e-5** | **7.7e-6** | 6.8e-5 | 1.0e-5 |
| exponential α = 8 | 3.9e-3 | 1.6e-4 | 4.8e-4 | 1.9e-6 |
| shifted exp (0 at the cutoff), α = 1 | 9.2e-4 | 7.9e-5 | 1.9e-4 | **1.8e-9** |
| Wendland C2 (from `e1`, N = 800 / 128²) | 2.3e-3 | 2.7e-4 | 4.5e-4 | 1.3e-5 |

- The exponential weight is 0.37 at the cutoff for α = 1; in 2-D, jitter moves neighbours across it,
  which adds a jump to the defects (superlinear growth with rf). A weight that vanishes at the cutoff
  removes that part, but the bulk of the defect is the general MLS asymmetry, so the weight shape alone
  is not a cure; Wendland is worse.
- α ≈ 2–4 with the library exponential reduces the 1-D drift 20×; α = 8 is too narrow.
- Larger R reduces the drift (and the defects) roughly as 1/R; R = 5.5 with the shifted exponential is
  almost exactly conservative in this 2-D test.

## 3. Changing the scheme

All variants in the standalone solver; MUSCL2 unless stated (`e2_forms.jl`, `e7_intercept.jl`,
`e8_local.jl`, `e3_ho_smooth.jl`).

| form | idea | conservation | consistency | cost |
|---|---|---|---|---|
| `gfd` (library) | −2Σ c_ij·(F*_ij − F_i) | defects above | 2nd order | – |
| `gfd` + intercept fit | MLS **with** a constant term over the neighbours (not the particle) ⇒ Σ_j c_ij = 0 exactly | row-sum defect removed; 1-D drift 3–15× smaller (1e-3 → 6.5e-5…1.9e-4), 2-D mixed | unchanged (2nd order, same errors) | none (local) |
| `sym` | face form −(1/V_i)Σ A_ij·F*_ij, A_ij = V_i c_ij − V_j c_ji | **exact** | **inconsistent**: smooth error plateaus at 1.4e-2…2.7e-2 (rf = 0.2) | local |
| `symc` | `sym` + antisymmetric correction with Σ_j A_ij = 0 (graph-Laplacian solve) | exact | zeroth order only: plateau ~1e-3 | one sparse solve |
| `symc1` | minimal-norm antisymmetric correction with Σ_j A_ij = 0 **and** Σ_j A_ij ⊗ d_ij = 2V_i I | **exact** | **2nd order**, smooth errors as `gfd` at R = 5.5 (1-D, 6.7e-5 vs 5.7e-5 at N = 800), 2–3× larger at R = 2.5; post-shock errors equal (1-D, N = 800, BJ: 1.2e-3 vs 1.1e-3 at R = 2.5, 2.2e-3 both at R = 5.5; 2-D within 10–25 %) | one sparse solve per geometry (2-D: 0.5 s for 4k, 2.9 s for 16k, 28 s for 65k particles, direct LU) |

Notes on `symc1`:
- The correction can flip the sign of A_ij relative to x_j − x_i; with the library's per-axis upwinding
  that is anti-dissipative and blows up (1-D, R = 2.5). Solving the Riemann problem along the face
  normal A_ij/|A_ij| fixes it: `A_ij·F* = |A_ij| f̂(n; u_i, u_j)` is dissipative for any A.
- With Σ_j A_ij = 0 and the normal Riemann solver, the first-order scheme is a convex combination for
  any coefficients: **maximum principle / positivity of the first-order fallback**, which the fixed-MLS
  GFD form does not have (`e9_maxprinciple.jl`: first-order Burgers box, min/max of u ∈ [0,1]: GFD with
  the cubic divergence −5.4e-2 / 1.004; face form 0 / 0.996, mass exact). This is the clainshock4 problem
  (MOOD fallback with MLS_order = 4) in scalar form: a conservative face form with a fixed geometry
  would give conservation and a positive fallback at the same time.
- Higher order: imposing the moment conditions up to degree p (`face_coefficients_p`) gives ~3rd order
  for p = 3 (MUSCL4 reaches 4 in the `gfd` form); a high-order conservative variant needs more work.

## 4. Side results

**Why MUSCL3 and MUSCL5 lose an order** (`e4_oddorder.jl`, `e5_divdegree.jl`). Not the Rusanov
dissipation (switching it off changes nothing). It is the divergence fit: the derivative of a
degree-q least-squares fit on a (near-)symmetric stencil is O(h^q) for even q (the next monomial d^{q+1}
is odd and aliases into the gradient) and O(h^{q+1}) for odd q (the next monomial is even and
orthogonal). MUSCL(o) uses q = o − 1, so even orders gain one, odd orders do not. Decoupling confirms it
(smooth advection, regular and 20 % jittered):

| reconstruction degree | divergence degree | EOC |
|---|---|---|
| 2 (MUSCL3) | 2 (now) | 2.0 |
| 2 (MUSCL3) | **3** | **3.0** |
| 4 (MUSCL5) | 4 (now) | ~4 |
| 4 (MUSCL5) | **5** | **~5** (with the lower regularisation below) |
| 3 | 2 | 2.0 |

So MUSCL(odd) needs a divergence fit of degree o (one above the reconstruction); the library caps
`div_order` at the scheme order, which would need to change. On 20 % jitter the odd-degree bonus fades
(MUSCL4: EOC 3.6 at N = 400).

**The ~1e-8 error floor of the high orders** is the moment-matrix regularisation (1e-10·tr(M)/B):
with 1e-14 the degree-5 divergence converges to 2.6e-11 (1-D, N = 800). Float32 needs a larger value,
so it should scale with eps(T).

## 5. Recommendations

1. Merge the seam fix (HyperCloudCore PR #9).
2. Keep the fixed MLS order with MOOD (coupling costs ~100× in mass).
3. Cheap: intercept MLS fit for the divergence (local, Σ_j c_ij = 0 exactly; 3–15× less drift in 1-D,
   no accuracy loss), exponential α ≈ 2–4 or a weight that vanishes at the cutoff, R ≥ 3.5.
4. Exact: the `symc1` face form (antisymmetric coefficients with zero row sums and linear exactness,
   normal-direction Riemann solver): exactly conservative, second order, positive first-order fallback,
   one sparse solve per geometry (static grids). Moving grids would need an iterative, warm-started solve.
5. Odd orders: divergence degree o for MUSCL(o) with odd o. Regularisation ∝ eps(T).

## 6. Implemented (HyperCloudCore PR #10, stacked on #9)

Recommendations 3 (intercept fit) and 5 (odd orders, regularisation) are in the library. The odd-order
bump only applies where the stencil has ≥ 1.5× as many neighbours as unknowns (1-D MUSCL3 needs
R ≥ 3.5Δx: an exactly determined degree-3 fit is fragile on moving grids), up to degree 3 in 3-D.
Library before (`env_base`, seam fix only) vs after (`env_study`), `e10_library_after.jl`, `out/e10.txt`
(periodic, 20 % jitter):

| | before | after |
|---|---|---|
| 1-D MUSCL1, R 2.5 / 5.5, drift after the shock | 6.6e-5 / 1.9e-5 | 4.7e-6 / 1.1e-6 |
| 1-D MUSCL2 VK, R 2.5 / 5.5 | 1.0e-3 / 3.4e-4 | 1.1e-4 / 7.5e-5 |
| 2-D MUSCL2 VK, R 2.5 / 5.5 | 1.7e-4 / 1.7e-6 | 7.0e-5 / 7.9e-6 |
| 1-D MUSCL4 MOOD MLS 4, R 2.5 / 5.5 | 9.4e-3 / 2.2e-4 | 1.8e-2 / 3.4e-5 |
| 1-D MUSCL3 MOOD MLS 3, R 2.5 / 5.5 | 4.9e-4 / 9.8e-5 | 3.7e-4 / 3.8e-4 |
| 2-D MUSCL3 MOOD MLS 3, R 2.5 / 5.5 | 2.8e-5 / 3.5e-5 | 2.3e-5 / 3.8e-6 |
| smooth MUSCL3 EOC, 1-D R 3.5 / 2-D R 2.5 | 2.0 / 2.1 | 3.0 / 3.0 |
| smooth MUSCL4 (2-D R 2.5), L1 at 128² | 4.2e-7 | 2.4e-7 |

Worse: 1-D MUSCL4 at R = 2.5 (4 neighbours for 4 unknowns: an exactly determined fit; that setup was
already poor, and also fails on moving grids before the change), 1-D MUSCL3 MOOD at R = 5.5 (degree-3
fits are more asymmetric than degree 2), 2-D MUSCL2 at R = 5.5 (tiny either way). The 2-D Euler
quadrant (400²) and DMR (960×300) runs on the GPU (Float32) are unchanged in extrema and step counts.
