# Isentropic vortex for the 2-D discrete-velocity BGK model (γ = 1.4: Chu-reduced pair, 2 Nv² equations), periodic
# [0, 10]², background (ρ, u, v, p) = (1, 1, 1, 1), one period. Near the Euler limit (small τ) the exact solution after
# t = 10 is the initial condition. Particle-management study 2026-10-28 (tools/adapt_study.jl, ADAPT_VORTEX_ERROR).
isdefined(@__MODULE__, :vortex_state) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "Euler", "IsentropicVortex.jl"))

struct BGKVortexIC{T} <: HyperCloud.InitialCondition end
(::BGKVortexIC{T})(pos) where {T} = State{4, T}(T.(vortex_state(Float64(pos[1]), Float64(pos[2]))))
HyperCloud.build_ic(::Val{:bgk_vortex}, conf::Dict, ctx::Dict) = BGKVortexIC{ctx[:T]}()

function get_config()
    shared = Dict{Symbol, Any}(
        :PDE_name => :bgk, :PDE_D => 2, :PDE_gamma => 1.4,
        :PDE_Nv => 8, :PDE_vmax => 7.0, :PDE_velocity_grid => :gauss_hermite,
        :PDE_tau => 1e-3,
        :sim_func_name => :run_bgk_simulation,
        :T => Float64, :tmax => 10.0, :snapshots => 10,
        :Grid_domain => :rectangular, :Grid_mins => (0.0, 0.0), :Grid_maxs => (10.0, 10.0), :Grid_Ns => (50, 50),
        :Grid_periodic => true, :Grid_randomness_factor => (0.2, 0.2), :Grid_seed => 42, :Grid_mover => :none,
        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 2.5,
        :IC_name => :bgk_vortex,
        :Time_CFL => 0.4,
    )
    methods = Dict{Symbol, Dict{Symbol, Any}}(
        :ARS222MUSCL2 => Dict{Symbol, Any}(:Time_stepper => :ARS222, :Scheme_name => :MUSCL, :Scheme_order => 2,
            :Scheme_MLS_order => 2, :Flux_name => :Rusanov, :Limiter_name => :none, :MOOD_criterion => :none),
        :ARS233MUSCL3 => Dict{Symbol, Any}(:Time_stepper => :ARS233, :Scheme_name => :MUSCL, :Scheme_order => 3,
            :Scheme_MLS_order => 2, :Flux_name => :Rusanov, :Limiter_name => :none, :MOOD_criterion => :none),
    )
    return SimulationConfig(shared, methods, collect(keys(methods)); sources = [@__FILE__])
end
