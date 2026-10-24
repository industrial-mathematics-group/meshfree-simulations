# Mass-conservation study data (Float64, CPU runs only, 2026-10-24)

Flat copy (`data/mass/data/*.jld2`, save-path layout) of the runs of the mass-conservation study, made with
`cluster/run_mass.sh` (one method, one parameter set, one N range per job; `SIM_SET` overrides the jitter,
weight range, alpha and MLS order of the experiment files, so these runs are **not** the tasks of the experiment
files themselves). Statistics from the hyperbolic stat preset (`mass`, `relative_mass`, `mass_error`, …).

| Study | Experiment | Parameters | N | Jobs |
|---|---|---|---|---|
| S1 orders × limiting | conservation_burgers (1-D shock, outflow) | jitter 0.4, R 5.5, α 1; orders 2–5 × {plain, EPD1, EPD2, VK, minmod}, fixed MLS = order, and EPD1/EPD2 with MLS 0 | 32–3163 | 24610696–24610729; minmod order 2 N = 1996, 2512, 3163: 24612208 24612209 24612210 |
| S2 jitter | conservation_burgers | jitter 0, 0.1, 0.2, 0.3, 0.5; orders 2, 4; plain and EPD1 | 32–3163 | 24611115–24611272 |
| S3 range, alpha | conservation_burgers | R 3.5, 4.5, 7.5, 9.5 (orders 2, 4, EPD1); α 0.5, 2, 4 (order 4, EPD1) | 32–3163 | 24611276–24611287 |
| S4 smooth periodic | Convergence1D (1-D linear advection) | jitter 0.2, 0.4; orders 2–5 plain | 32–3163 | 24611289–24611296 |
| S5 2-D | ConvergenceBurgers (2-D periodic, pre-shock) | jitter 0, 0.2, 0.4 × R 3.5, 5.5 × orders 2, 4 × {plain, EPD1, VK} | 32–200 | 24611297–24611342 |

The exact job of a run: `data/<jobid>/job-info.txt` next to its copy. Tables: `reports/2026-10-24-mass-conservation.md`
(summary and interpretation) and `reports/2026-10-24-mass-conservation-tables.md` (every run), both from
`julia --project=env_headless tools/mass_table.jl data/mass`; CSV of all runs in `mass_runs.csv` here.
