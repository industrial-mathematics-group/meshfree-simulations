# Noh problem in 2-D (Noh 1987; Liska & Wendroff 2003): cold gas (ρ = 1, p = p0 = 1e-6) streams radially inwards with
# unit speed onto the origin, γ = 5/3, [-1, 1]², t = 0.6. Exact solution: a cylindrical shock at r = t/3, behind it
# ρ = 16, u = 0, p = 16/3; in front ρ = 1 + t/r, u = −x/r, p = p0 (1 + t/r)^γ. Strong shock, and "wall heating" at
# the centre (density dip / temperature spike of most schemes).
# The exact state is prescribed in the ghost layer at every time (NohExactBC) and is the reference of the error
# statistics of the hyperbolic preset (`reference = :analytical_solution` with the `analytic_closure` below).
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

"Exact Noh state (ρ, ρu, ρv, E) at (x, y) and time t."
@inline function noh_exact(x::T, y::T, t::T, γ::T, p0::T) where {T}
    r = sqrt(x^2 + y^2)
    if r < (γ - 1) / 2 * t                       # behind the shock (speed (γ − 1)/2 = 1/3 for γ = 5/3), gas at rest:
        ρ = ((γ + 1) / (γ - 1))^2                # ρ = 16
        p = (γ - 1) * ρ / 2                      # internal energy 1/2 (the kinetic energy of the inflow): p = 16/3
        return State{4, T}(ρ, zero(T), zero(T), p / (γ - 1))
    end
    ρ = r > 0 ? 1 + t / r : one(T)
    ux = r > 0 ? -x / r : zero(T); uy = r > 0 ? -y / r : zero(T)
    p = p0 * ρ^γ
    return State{4, T}(ρ, ρ * ux, ρ * uy, p / (γ - 1) + ρ * (ux^2 + uy^2) / 2)
end

struct NohIC{T} <: HyperCloud.InitialCondition
    γ::T
    p0::T
end
(ic::NohIC{T})(pos) where {T} = noh_exact(T(pos[1]), T(pos[2]), zero(T), ic.γ, ic.p0)
HyperCloud.build_ic(::Val{:noh}, conf::Dict, ctx::Dict) = NohIC{ctx[:T]}(ctx[:T](conf[:gamma]), ctx[:T](conf[:p0]))

HyperCloud.analytic_closure(eq::HyperCloud.EulerEquation{2, 4, T}, ic::NohIC, geom::GeometricDomain) where {T} =
    st -> noh_exact(T(st[1]), T(st[2]), T(st[end]), T(eq.gamma), T(ic.p0))

"Ghost states from the exact solution at time t (all boundaries)."
struct NohExactBC <: AbstractBoundaryCondition
    γ::Float64
    p0::Float64
end
@inline function _noh_bc!(i, buf, pos, tags, is_boundary, tag, t, γ, p0)
    if is_boundary[i] && tags[i] == tag
        buf[i] = noh_exact(pos[i][1], pos[i][2], t, γ, p0)
    end
    return nothing
end
function (bc::NohExactBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    T = eltype(eltype(rhos_buffer))
    foreach_particle(_noh_bc!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.positions, pg.core.tags, pg.core.is_boundary, tag, T(t), T(bc.γ), T(bc.p0)))
    return nothing
end

function get_config()
    γ, p0 = 5 / 3, 1e-6
    bc = NohExactBC(γ, p0)
    shared = bench_shared(:PDE_gamma => γ, :tmax => 0.6, :Grid_mins => (-1.0, -1.0), :Grid_maxs => (1.0, 1.0),
                          :Grid_Ns => (400, 400), :Grid_bc_map => Dict{Int, AbstractBoundaryCondition}(k => bc for k in 1:4),
                          :IC_name => :noh, :IC_gamma => γ, :IC_p0 => p0)
    return bench_config(shared; reference = :analytical_solution, sources = [@__FILE__])
end
