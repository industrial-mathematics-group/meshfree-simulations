# Convergence data (Float64, CPU runs only)

Flat copy of the *correct* convergence runs of 2026-10-22/23 (`data/convergence/data/*.jld2`, the save-path
layout: `set_save_path!("data/convergence")` or point `SIM_SAVE` at this folder and the runs are found as cached).
All runs: Float64, CFL 0.2, R = 5.5, regular lattice (jitter 0), RK4, analytical reference, HyperCloud 0.5.0 /
HyperCloudCore 0.7.0 / PDEStudioCore 0.4.1. Everything here ran on the cluster's CPU nodes; the V100 runs of the
same sweeps (jobs 24544694–697, 24544709–711, 24544716, 24558179/180) are **not** included because the CUDA path
computes partly in single precision (error floor ~2e-7 mean for orders ≥ 4, see
reports/2026-10-22-convergence-clainshock.md §1b). If you have their files in your own data folder, they are the
ones to replace: 2-D advection and Burgers at N ≥ 251 (252, 317, 399, 502, 631, 795, 1000).

| Experiment | N | Jobs (CPU) |
|---|---|---|
| Convergence1D, 4 methods | 32 … 3163 (21) | 24544690–693 |
| ConvergenceAdvection, 4 methods | 32 … 200 (9) | 24544702–705 |
| ConvergenceAdvection MUSCL4 / MUSCL5 | 252 … 1000 (7) | 24559671 / 24559672 + 24559673 |
| ConvergenceAdvection MUSCL2 / MUSCL3 | 252 … 1000 (7) | 24588423 / 24588424 |
| ConvergenceBurgers, 4 methods | 32 … 200 (9) | 24544823, 24544831, 24544844, 24544857 |
| ConvergenceBurgers MUSCL2 / MUSCL3 / MUSCL4 / MUSCL5 | 252 … 1000 (7) | 24588466 / 24588467 / 24588468 / 24588469 + 24588470 |

The table of errors and orders: `reports/2026-10-23-convergence-table.md` (made by
`julia --project=env_headless tools/convergence_table.jl data/convergence`), CSV in `convergence_table.csv` here.
