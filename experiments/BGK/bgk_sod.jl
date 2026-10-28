# BGK Sod shock tube (HyperCloud's bgk_1d_sod) with the methods of the particle-management study added. Moments are
# compared with the exact Euler solution for small τ.
function get_config()
    config = bgk_1d_sod()
    config.shared[:T] = Float64
    for (st, K) in ((:ARS222, 2), (:ARS233, 3))
        config.methods[Symbol("$(st)MUSCL$(K)(EPD2)")] = Dict{Symbol, Any}(:Time_stepper => st, :Scheme_name => :MUSCL,
            :Scheme_order => K, :Scheme_MLS_order => 2, :Flux_name => :AlignedRusanov, :Limiter_name => :none,
            :MOOD_criterion => :U2, :MOOD_strategy => :EPD2, :MOOD_delta_relax => 1e-4)
    end
    push!(config.sources, @__FILE__)
    return config
end
