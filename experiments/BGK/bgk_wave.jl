# BGK density wave: ρ = 1 + 0.2 sin(2πx), u = 1, p = 1 on the periodic unit interval (γ = 1.4: Chu-reduced pair on a
# 1-D velocity grid). In the Euler limit (τ → 0) the exact solution is the translated profile; for any τ the
# Maxwellian with constant u and p is advected without change in the Euler limit, and a grid moving with u = 1
# (Lagrangian) carries it without relative motion. Particle-management study 2026-10-28 (tools/adapt_study.jl).
const WAVE_U = 1.0
wave_density(x) = 1 + 0.2 * sinpi(2x)

struct WaveIC{T} <: HyperCloud.InitialCondition end
function (::WaveIC{T})(pos) where {T}
    ρ = wave_density(Float64(pos[1])); u = WAVE_U; p = 1.0
    State{3, T}(T.((ρ, ρ * u, p / 0.4 + ρ * u^2 / 2)))
end
HyperCloud.build_ic(::Val{:density_wave}, conf::Dict, ctx::Dict) = WaveIC{ctx[:T]}()

function get_config()
    shared = Dict{Symbol, Any}(
        :PDE_name => :bgk, :PDE_D => 1, :PDE_gamma => 1.4,
        :PDE_Nv => 24, :PDE_vmax => 7.0, :PDE_velocity_grid => :uniform,
        :PDE_tau => 1e-4,
        :sim_func_name => :run_bgk_simulation,
        :T => Float64,
        :tmax => 1.0, :snapshots => 10,
        :Grid_domain => :rectangular, :Grid_mins => (0.0,), :Grid_maxs => (1.0,), :Grid_Ns => (200,),
        :Grid_periodic => true, :Grid_randomness_factor => (0.2,), :Grid_seed => 42, :Grid_mover => :none,
        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 5.5,
        :IC_name => :density_wave,
        :Time_CFL => 0.2,
    )
    methods = Dict{Symbol, Dict{Symbol, Any}}()
    for (st, K) in ((:ARS222, 2), (:ARS233, 3))
        methods[Symbol("$(st)MUSCL$(K)")] = Dict{Symbol, Any}(:Time_stepper => st, :Scheme_name => :MUSCL, :Scheme_order => K,
            :Scheme_MLS_order => 2, :Flux_name => :Rusanov, :Limiter_name => :none, :MOOD_criterion => :none)
        methods[Symbol("$(st)MUSCL$(K)(VK)")] = Dict{Symbol, Any}(:Time_stepper => st, :Scheme_name => :MUSCL, :Scheme_order => K,
            :Scheme_MLS_order => 2, :Flux_name => :Rusanov, :Limiter_name => :VK, :Limiter_mode => :hard, :MOOD_criterion => :none)
        methods[Symbol("$(st)MUSCL$(K)(EPD2)")] = Dict{Symbol, Any}(:Time_stepper => st, :Scheme_name => :MUSCL, :Scheme_order => K,
            :Scheme_MLS_order => 2, :Flux_name => :AlignedRusanov, :Limiter_name => :none,
            :MOOD_criterion => :U2, :MOOD_strategy => :EPD2, :MOOD_delta_relax => 1e-4)
    end
    return SimulationConfig(shared, methods, collect(keys(methods)); sources = [@__FILE__])
end
