# Simulations

## Which config belongs to which data

All runs of the current data set (HyperCloudCore 0.7.8 rerun, `reports/2026-10-29-rerun-078.md`) are stored flat in
`data/*.jld2`, the save path of `startup.jl local` (the repository root). The configs below reproduce the exact
parameters of these runs (the job scripts in `cluster/` override precision, grid, CFL, sweep and methods), so
`@sim <config>` loads every run from `data/` instead of recomputing it.

| Config (`@sim …`) | Runs | Backend, precision |
|---|---|---|
| `Mach3Step_900x300` | 19 methods, 900×300, CFL 0.8 | GPU, Float32 |
| `DMR_1440x450`, `DMR_960x300` | 17 methods / the 2 coupled K = 4 methods | GPU, Float32 |
| `clainshock_800x800`, `clainshock_outflow_800x800` | 19 methods (exact boundaries) / 3 methods (outflow boundaries) | GPU, Float32 |
| `IsentropicVortex_Float32` | 6 methods, N = 50 … 400 | GPU, Float32 |
| `IsentropicVortex_Float64` | RK4SSPMUSCL4(none,AR,MLS4), N = 50 … 400 (precision check of the Float32 floor) | CPU (local), Float64 |
| `Convergence1D_sweep`, `ConvergenceAdvection_sweep`, `ConvergenceBurgers_sweep` | K = 2 … 5, lattice and jitter 0.2, R = 5.5 / 2.5 | CPU, Float64 |
| `defect_floor`, `defect_shock`, `defect_jitter`, `defect_range`, `defect_alpha`, `defect_smooth1d`, `defect_smooth2d` | conservation-defect experiments | CPU, Float64 |
| `bgk_sod_sweep` | 3 methods × static / Lagrangian / Lagrangian + adaptivity × τ × N (one run failed, see the report) | CPU, Float64 |

- `data/vortex-cpu-float32/` holds the same vortex runs computed on the CPU in Float32 (same parameters, hence the
  same hashes as the GPU runs): `set_save_path!("data/vortex-cpu-float32")` to load those instead.
- `data/jobinfo/` has the cluster job info (job id, node, GPU) of every fetched job.
- The other files in `experiments/` are the base configs the ones above build on (`Mach3Step.jl`, `DMR.jl`,
  `clainshock.jl`, `clainshock_outflow.jl`, `IsentropicVortex.jl`, `Convergence*.jl`, `conservation_burgers.jl`,
  `bgk_sod.jl`) and `methods_review.jl` (the method set, labels `<stepper>MUSCL<K>(<variant>,<AR|RU>,MLS<m>)`).
  `@sim` of a base config uses its own defaults and computes new runs.
- `bgk_wave.jl`, `bgk_vortex.jl` and the particle-management cases run through `tools/adapt_study.jl` (own time loop
  with per-call diagnostics); their results are the RESULT lines in `reports/`, not `@sim` runs.

Check that a config finds all its runs (and which parameters differ for files it does not match):

    julia --project=env_headless tools/check_config.jl [--orphans] Mach3Step_900x300 DMR_1440x450 ...
