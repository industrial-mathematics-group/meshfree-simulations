# MOOD u2 diagnostics (2026-10-28)

Scripts behind `reports/2026-10-28-mood-u2-vortex.md`. They run with an environment that has the registered
HyperCloud and a HyperCloudCore checkout (`Pkg.develop(path = …)`), created outside the repo; run them from this
folder (`include("u2_override.jl")` is relative).

| Script | What it does |
|---|---|
| `vortex_mood.jl N K [tmax] [count]` | isentropic vortex (review setup), RK/MUSCL K with EPD2/U2; with `count = 1` it replaces the u2 criterion by an instrumented copy (`u2_override.jl`) that counts every DMP violation and why u2 rejected it, and histograms position, overshoot and local range. Env: `U2_VARIANT` (letters a = units, b = violating components only, c = no ratio test, p = plateau), `U2_ETA`, `DELTA_RELAX`, `VORTEX_NOMOOD`. Use it on HyperCloudCore **main** (≤ 0.7.5): it overrides `detect_and_drop!`. |
| `shock_mood.jl CASE K N` | discontinuity checks with the same override: `sod` (1-D, N particles) or `box` (2-D advection, N²); prints L1 and over/undershoot. `DELTA_ABS` sets an absolute δ. |
| `vortex_lib.jl`, `shock_lib.jl` | the same runs through the library criterion (HyperCloudCore `claude/mood-u2-extrema`): `LIB_PLATEAU`, `LIB_RATIO` set `:MOOD_plateau`, `:MOOD_curvature_ratio`. |
| `box1d.jl K` | HyperCloudCore's pipeline boundedness test (1-D box on a lattice) for several plateau/ratio settings. |
