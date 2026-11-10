# Opens the BGK Sod sweep (experiments/BGK/bgk_sod_sweep.jl, HyperCloudCore 0.8.0 runs) in a headless PDEStudio.
# Included by the experiment scripts of this folder. Environment: env_web (PDEStudio 0.6.1, HyperCloud 0.6.1).
# Data: data/pdestudio05/ (the 243 runs of 2026-11-09, upgraded to the PDEStudioCore 0.5 format).
using PDEStudio, HyperCloud, StaticArrays
const ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const IMG = joinpath(@__DIR__, "img"); mkpath(IMG)
set_target_module!(Main)
HyperCloud.CONFIG_ROOT_OVERWRITE[] = ROOT
set_save_path!(joinpath(ROOT, "data", "pdestudio05"))
set_stat_preset!("hyperbolic")
config = let c = Base.include(Main, HyperCloud.resolve_config_path(Main, "bgk_sod_sweep"))
    c isa SimulationConfig ? c : Base.invokelatest(Main.get_config)
end
const ALL_METHODS = copy(config.defaults)   # the studio edits config.defaults in place when methods change
studio_headless!()
launch_plotter(config)
shot(name; kw...) = studio_snapshot(joinpath(IMG, name * ".png"); kw...)
