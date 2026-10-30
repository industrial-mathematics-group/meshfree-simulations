# Shock–vortex interaction (Jiang & Shu 1996; Rault, Chiavassa & Donat 2003): a stationary Mach 1.1 shock at x = 0.5 in
# [0, 2] × [0, 1], γ = 1.4, upstream (ρ, u, v, p) = (1, 1.1 √γ, 0, 1), downstream from the Rankine–Hugoniot relations;
# an isentropic vortex centred at (0.25, 0.5) is superposed on the upstream flow and is convected through the shock
# (t = 0.6):  δu = ε τ e^{α(1−τ²)} sin θ,  δv = −ε τ e^{α(1−τ²)} cos θ,  δT = −(γ−1) ε² e^{2α(1−τ²)} / (4αγ),
# τ = r / r_c, ε = 0.3, r_c = 0.05, α = 0.204, with ρ = T^{1/(γ−1)}, p = ρ^γ in the vortex (T = p/ρ = 1 upstream).
# Tests the coupling of shock capturing with smooth vortical flow (the vortex is split, the shock is deformed).
# Boundaries: fixed upstream state on the left, outflow on the right, slip walls at the top and bottom.
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

struct ShockVortexIC{T} <: HyperCloud.InitialCondition
    γ::T
    uL::T            # upstream velocity
    R::NTuple{4, T}  # downstream conserved state
end
function (ic::ShockVortexIC{T})(pos) where {T}
    x, y = T(pos[1]), T(pos[2]); γ = ic.γ
    x >= T(0.5) && return State{4, T}(ic.R)
    ε, rc, α = T(0.3), T(0.05), T(0.204)
    dx, dy = x - T(0.25), y - T(0.5)
    τ = sqrt(dx^2 + dy^2) / rc
    e = exp(α * (1 - τ^2))
    sinθ, cosθ = τ > 0 ? (dy / (τ * rc), dx / (τ * rc)) : (zero(T), zero(T))
    u = ic.uL + ε * τ * e * sinθ
    v = -ε * τ * e * cosθ
    Tm = 1 - (γ - 1) * ε^2 * e^2 / (4 * α * γ)
    ρ = Tm^(1 / (γ - 1)); p = ρ^γ
    return State{4, T}(ρ, ρ * u, ρ * v, p / (γ - 1) + ρ * (u^2 + v^2) / 2)
end
function HyperCloud.build_ic(::Val{:shock_vortex}, conf::Dict, ctx::Dict)
    T = ctx[:T]; γ = conf[:gamma]; M = conf[:mach]
    uL = M * sqrt(γ)                                         # ρ = p = 1 upstream: sound speed √γ
    ρR = (γ + 1) * M^2 / ((γ - 1) * M^2 + 2)
    pR = (2γ * M^2 - (γ - 1)) / (γ + 1)
    uR = uL / ρR
    return ShockVortexIC{T}(T(γ), T(uL), T.(bench_cons(γ, ρR, uR, 0, pR)))
end

function get_config()
    bc = Dict{Int, AbstractBoundaryCondition}(1 => FixedDirichlet(), 2 => OutflowBC(), 3 => BenchWallY(), 4 => BenchWallY())
    shared = bench_shared(:tmax => 0.6, :Grid_mins => (0.0, 0.0), :Grid_maxs => (2.0, 1.0), :Grid_Ns => (800, 400),
                          :Grid_bc_map => bc, :IC_name => :shock_vortex, :IC_gamma => 1.4, :IC_mach => 1.1)
    return bench_config(shared; sources = [@__FILE__])
end
