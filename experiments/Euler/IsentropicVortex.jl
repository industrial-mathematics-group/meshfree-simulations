# Isentropic vortex (Yee, Sandham, Djomehri 1999): smooth Euler solution advected with the background velocity
# (1, 1) through the periodic box [0, 10]², so after t = 10 the exact solution is the initial condition again.
# Paper review 2026-10-26: shows whether MOOD leaves smooth flow at full order (run with MOOD on, also in single
# precision). The error is measured against the initial condition at t = tmax (tools/vortex_error.jl).
#   julia --project=env_headless tools/run_sim.jl IsentropicVortex     (SIM_VARY="Grid_Ns__all=[50, 100, 200]")
isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))

const VORTEX_GAMMA = 1.4
const VORTEX_BETA = 5.0
const VORTEX_CENTER = (5.0, 5.0)

"Conserved state (ρ, ρu, ρv, E) of the vortex at `x` (Float64), background (ρ, u, v, p) = (1, 1, 1, 1)."
function vortex_state(x, y)
    γ, β = VORTEX_GAMMA, VORTEX_BETA
    dx, dy = x - VORTEX_CENTER[1], y - VORTEX_CENTER[2]
    r2 = dx^2 + dy^2
    ex = exp((1 - r2) / 2)
    u = 1 - β / (2π) * ex * dy
    v = 1 + β / (2π) * ex * dx
    temp = 1 - (γ - 1) * β^2 / (8γ * π^2) * exp(1 - r2)
    ρ = temp^(1 / (γ - 1))
    p = ρ^γ
    E = p / (γ - 1) + ρ * (u^2 + v^2) / 2
    return (ρ, ρ * u, ρ * v, E)
end

struct VortexIC{T} <: HyperCloud.InitialCondition end
(::VortexIC{T})(pos) where {T} = State{4, T}(T.(vortex_state(Float64(pos[1]), Float64(pos[2]))))
HyperCloud.build_ic(::Val{:my_ic}, conf::Dict, ctx::Dict) = VortexIC{ctx[:T]}()

function get_config()
    varied_params = create_varied_dict()
    push!(varied_params, :Grid_Ns__all => [50, 71, 100, 141, 200, 283, 400])

    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 2,
        :PDE_gamma => VORTEX_GAMMA,
        :T => Float64,

        :sim_func_name => :run_direct_simulation,
        :tmax => 10.0,              # one period: the vortex returns to its initial position
        :snapshots => 10,

        :Grid_domain => :rectangular,
        :Grid_mins => (0.0, 0.0),
        :Grid_maxs => (10.0, 10.0),
        :Grid_Ns => (100, 100),
        :Grid_periodic => true,
        :Grid_randomness_factor => (0.2, 0.2),
        :Grid_seed => 42,
        :Grid_mover => :none,

        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,       # the range of the shock problems (the review's point 3)

        :IC_name => :my_ic,

        :Time_CFL => 0.8,
        :Scheme_MLS_order => 2,
    )

    methods = get_master_method_dict()
    add_review_methods!(methods)
    defaults = ["RK4SSPMUSCL4(EPD2,AR,MLS2)", "RK4SSPMUSCL4(EPD2,AR,MLS4)", "RK4SSPMUSCL4(none,AR,MLS4)",
                "RK3MUSCL3(EPD2,AR,MLS2)", "RK3MUSCL3(none,AR,MLS3)", "RK2MUSCL2(EPD2,AR,MLS2)"]
    return SimulationConfig(shared, methods, defaults; varied = varied_params, sources = [@__FILE__])
end
