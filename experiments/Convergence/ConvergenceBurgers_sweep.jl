# ConvergenceBurgers_sweep: the convergence sweeps of the HyperCloudCore 0.7.8 rerun (reports/2026-10-29-rerun-078.md §1), CPU, Float64:
# SSP stepper, AlignedRusanov flux, no limiting, MLS order K, lattice and jitter 0.2, R = 5.5 (K = 2 … 5) and
# R = 2.5 (K = 2 … 4), α = 1. Each (K, jitter, R) is a method copy with the parameters cluster/run_mass.sh set, so
# `@sim ConvergenceBurgers_sweep` finds every run as cached.
# Check: julia --project=env_headless tools/check_config.jl ConvergenceBurgers_sweep
module _ConvergenceBurgers_sweep_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "ConvergenceBurgers.jl"))
end

function get_config()
    config = _ConvergenceBurgers_sweep_Base.get_config()
    config.shared[:T] = Float64
    empty!(config.varied)
    config.varied[:Grid_Ns__all] = collect(ceil.(Int, 10 .^ (1.5:0.1:3.0)))
    config.defaults = Symbol[]
    for K in 2:5, jit in (0.0, 0.2), R in (5.5, 2.5)
        R == 2.5 && K == 5 && continue
        base = Symbol(_ConvergenceBurgers_sweep_Base.review_label(K, :none, :AlignedRusanov, K))
        name = Symbol("$(base) jitter$(jit) R$(R)")
        config.methods[name] = merge(config.methods[base], Dict{Symbol, Any}(
            :Grid_randomness_factor => (jit, jit), :Weight_range => R, :Weight_alpha => 1.0, :Scheme_MLS_order => K))
        push!(config.defaults, name)
    end
    push!(config.sources, @__FILE__)
    return config
end
