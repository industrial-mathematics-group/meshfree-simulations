# Convergence data on jittered points (randomness factor 0.2, Float64, CPU runs only)

Same three experiments and parameters as `data/convergence/` (regular lattice), but with
`Grid_randomness_factor` 0.2 in every direction (`cluster/run_sweep.sh ... jitter0.2`, i.e. the same jitter the
Euler cases use), run 2026-10-23/24 on the cluster's CPU nodes in Float64. Flat save-path layout:
`set_save_path!("data/convergence_jitter02")` finds the runs as cached. Note that these runs are *not* the
tasks of the experiment files (the configs have jitter 0): they carry `:Grid_randomness_factor => (0.2, 0.2)`
(`(0.2,)` in 1-D).

| Experiment | N | Jobs |
|---|---|---|
| Convergence1D, 4 methods | 32 … 3163 | 24596733 24596734 24596735 24596736 |
| ConvergenceAdvection, N ≤ 200 / N ≥ 251 | 4 methods | 24596737 24596738 24596739 24596740 24596741 24596742 24596743 24596744 24596756 |
| ConvergenceBurgers, N ≤ 200 / N ≥ 251 | 4 methods | 24596757 24596758 24596759 24596760 24596775 24596776 24596777 24596778 24596779 |

Table: `reports/2026-10-24-convergence-table-jitter02.md` (`tools/convergence_table.jl`), CSV here; comparison
with the lattice set: `python3 tools/compare_convergence.py data/convergence/convergence_table.csv
data/convergence_jitter02/convergence_table.csv lattice jitter0.2`.
