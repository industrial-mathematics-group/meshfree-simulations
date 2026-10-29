# bgk_sod_sweep: the BGK Sod runs of the HyperCloudCore 0.7.8 rerun (reports/2026-10-29-rerun-078.md §6), CPU, Float64:
# ARS222MUSCL2(VK), ARS222MUSCL2(EPD2), ARS233MUSCL3(EPD2) on a static grid, on a Lagrangian grid and on a Lagrangian
# grid with particle management (MLS volumes), τ = 1e-5, 1e-3, 1e-2, N = 100 … 1600. Each (method, grid) is a method
# copy with the parameters cluster/run_bgk.sh set, so `@sim bgk_sod_sweep` finds every run as cached.
# Check: julia --project=env_headless tools/check_config.jl bgk_sod_sweep
module _bgk_sod_sweep_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "bgk_sod.jl"))
end

function get_config()
    config = _bgk_sod_sweep_Base.get_config()
    config.shared[:T] = Float64
    empty!(config.varied)
    config.varied[:Grid_Ns__1] = [100, 141, 200, 283, 400, 566, 800, 1131, 1600]
    config.varied[:PDE_tau] = [1e-5, 1e-3, 1e-2]
    config.defaults = Symbol[]
    for base in ("ARS222MUSCL2(VK)", "ARS222MUSCL2(EPD2)", "ARS233MUSCL3(EPD2)"),
        (grid, mover, adapt) in (("static", :none, false), ("lagrangian", :lagrangian, false), ("lagrangian+adapt", :lagrangian, true))
        name = Symbol("$base $grid")
        config.methods[name] = merge(config.methods[Symbol(base)], Dict{Symbol, Any}(
            :Grid_mover => mover, :Grid_volumes => :mls, :Grid_volumes_degree => 2, :Adapt_enabled => adapt))
        push!(config.defaults, name)
    end
    push!(config.sources, @__FILE__)
    return config
end
