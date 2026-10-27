# ConvergenceAdvection_review: the convergence sweeps of the 2026-10-26 setup (reports/2026-10-27-review-rerun.md §1): SSP stepper,
# AlignedRusanov flux, no limiting, MLS order K, Float64, lattice and jitter 0.2, R = 5.5 and 2.5 (K <= 4 at 2.5).
# Each (K, jitter, R) is a method copy with these overrides, so the fetched results (data/review/convergence) are
# found as cached runs.
module _ConvergenceAdvection_review_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "Convergence", "ConvergenceAdvection.jl"))
end

function get_config()
    config = _ConvergenceAdvection_review_Base.get_config()
    config.shared[:T] = Float64
    config.defaults = Symbol[]
    for K in 2:5, jit in (0.0, 0.2), R in (5.5, 2.5)
        R == 2.5 && K == 5 && continue
        base = Symbol(_ConvergenceAdvection_review_Base.review_label(K, :none, :AlignedRusanov, K))
        name = Symbol("$(base) jitter$(jit) R$(R)")
        config.methods[name] = merge(config.methods[base], Dict{Symbol, Any}(
            :Grid_randomness_factor => (jit, jit), :Weight_range => R, :Weight_alpha => 1.0, :Scheme_MLS_order => K))
        push!(config.defaults, name)
    end
    push!(config.sources, @__FILE__)
    return config
end
