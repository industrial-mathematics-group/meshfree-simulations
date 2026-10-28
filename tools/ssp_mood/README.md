# MOOD substep checks and the K = 1 orientation (2026-10-28)

Scripts behind `reports/2026-10-28-ssp-mood-substeps.md`.

- `ssp_check.jl`: canonical Shu–Osher forms of the explicit tableaus (SSP coefficient, which substeps each stage
  uses); runs with any HyperCloudCore environment.
- `ssp_study.jl CASE K N SUBSTEPS [PLATEAU]`: one run (cases `sod`, `r123`, `vac`, `box`, `vortex`), printing L1,
  bounds/positivity and wall time; `:Time_mood_substeps` is set from SUBSTEPS (needs HyperCloudCore with the key;
  older versions ignore it). Env `CFL` for the Riemann cases.
- `ssp_batch.sh`, `vac_batch.sh`: the batches of the report. They `cd` into the folder they were written in and use
  two throwaway environments (`env_ssp`: HyperCloud + HyperCloudCore `claude/ssp-mood-substeps`; `env_main2`:
  HyperCloud + HyperCloudCore main); adjust both paths before rerunning.
