# 2026-11-10: BGK figures with PDEStudio 0.6's agent API

Figures of the BGK Sod runs of the HyperCloudCore 0.8.0 rerun (reports/2026-11-09-bgk-rerun-080.md), made headless
with PDEStudio 0.6.0's agent API (`studio_headless!`, `describe_studio`, `studio_set!`, `studio_option!`,
`studio_view!`, `studio_snapshot`, `studio_data`). Script: `tools/studio/bgk_sod_figures.jl`; figures (gitignored):
`figures/bgk080/`.

## Setup

- **Versions:** PDEStudio 0.6.0 needs PDEStudioCore 0.5, HyperCloud 0.6.0 still needs PDEStudioCore 0.4, so they can't
  share an environment. The figures use a scratch environment (PDEStudio 0.6.0, HyperCloudCore 0.8.0 for the boundary
  condition types, CairoMakie) outside the repo; the repo environments are unchanged.
- **Data migration:** the 243 BGK Sod runs of 2026-11-09 were copied to `data/pdestudio05/` and upgraded there with
  `upgrade_sim_data()` (243 upgraded, 0 unreadable). The originals in `data/` stay in the 0.4 format for the run
  environments. `upgrade_sim_data` does not rehash; lookups still match (below).
- **Config without HyperCloud:** `tools/studio/config_from_data.jl` rebuilds a `SimulationConfig` from the stored
  (normalised) parameters: constant ones shared, the swept keys varied, each other combination one method. All 243
  tasks find their file ("243 runs, 0 ran, 243 cached"), so nothing is recomputed and the solver is never needed.

## Figures

| File | Content |
|---|---|
| `sod_methods_N800_tau1e-5.png`, `…_zoom.png` | ρ at t = 0.2, N = 800, τ = 1e-5, all 9 variants (3 methods × static / Lagrangian / Lagrangian + particle management); zoom on contact and shock |
| `sod_tau_N800.png` | τ = 1e-5, 1e-3, 1e-2 side by side: VK static, ARS222 EPD2 Lagrangian, ARS233 EPD2 Lagrangian |
| `sod_convergence_l1.png` | final L1(ρ) against N per variant and τ (CairoMakie, see below) |

Observations:
- **ARS222MUSCL2(EPD2) on the Lagrangian grid loses the second plateau** (ρ ≈ 0.27 between contact and shock): the
  density ramps from x ≈ 0.15 to 0.34, with and without particle management. The same scheme on the static grid, VK on
  the Lagrangian grid and ARS233MUSCL3(EPD2) on the Lagrangian grid show both plateaus. Not new in 0.8.0 (L1 at
  N = 1600: 1.11e-2 on 0.7.8, 1.15e-2 now, vs 4.3e-3 static). Guess: second order + MOOD order drops on the clustered
  Lagrangian particles behind the contact.
- **τ:** at τ = 1e-2 the profile is smooth (far from equilibrium) and all methods agree; the errors against the Euler
  solution level off at ≈ 0.025 for every method (model error). At τ = 1e-5 static VK converges at first order (shock).

## PDEStudio 0.6 agent API: what worked, what didn't

Worked well: `describe_studio` lists every control with its valid values (an agent can act on it directly), invalid
names/values throw an `ArgumentError` with the choices, `studio_set!` applies the needed button press itself,
`studio_view!` zooms, `studio_snapshot` writes the image, `studio_state`/`studio_data` return the numbers. No display,
no browser, ~1 min per session including loading 243 runs.

Issues found:
1. **`particles_lines` with particle management:** lines join the particles in storage order; inserted particles are
   appended, so the lines fan across the plot. Lines should follow the position (sort per snapshot) for 1-D particles.
2. **Statistics only on the x axis:** in grid styles a statistic like `l1error|t` becomes the x coordinate while the
   dependent axis stays `Solution` (at the slider position): no error-over-time or error-vs-N plot. The U axis would
   need the statistics (`u = l1error` over `t` or over a swept parameter). In particle styles no statistic axes at all.
3. **`studio_set!` ordering:** `studio_set!(:compare_target => :none, :PDE_tau => 1e-5)` fails ("slider :PDE_tau is
   disabled"): the slider is checked before the layout change is applied. Two calls work.
4. **Small values in titles and sliders:** τ = 1e-5 shows as "0.000" (fixed 3 decimals); the swept key `Grid_Ns__1`
   is titled "Grid Ns  1", `PDE_tau` "Pde Tau".
5. **Long method names** overlap in titles at the default width (no wrapping or shrinking).

Verdict: **PASS** (all figures produced from migrated data without recomputation; findings above for PDEStudio).
