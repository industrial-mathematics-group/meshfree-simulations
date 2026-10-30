# Gresho vortex (Gresho & Chan 1990; Liska & Wendroff 2003; low-Mach variant as in Miczek, Röpke & Edelmann 2015):
# a steady rotating vortex, ρ = 1, γ = 1.4, centre (0.5, 0.5) of the periodic box [0, 1]², azimuthal velocity
#   u_φ = 5r (r < 0.2),  2 − 5r (0.2 ≤ r < 0.4),  0 (r ≥ 0.4),
# balanced by the pressure
#   p = p0 + 12.5 r²  (r < 0.2),  p0 + 4 − 4 ln 0.2 + 12.5 r² − 20 r + 4 ln r  (0.2 ≤ r < 0.4),  p0 − 2 + 4 ln 2  (r ≥ 0.4)
# (p is continuous). The initial condition is the exact solution at every time, so it is the reference of the error
# statistics. p0 sets the Mach number of the peak velocity 1: p0 = 5 is the classic setup (Ma ≈ 0.38),
# p0 = 1 / (γ Ma²) the low-Mach one, where compressible schemes without low-Mach correction are too dissipative.
# Configs: Gresho.jl, GreshoLowMach.jl.
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

@inline function gresho_state(x::T, y::T, γ::T, p0::T) where {T}
    dx, dy = x - T(0.5), y - T(0.5)
    r = sqrt(dx^2 + dy^2)
    if r < T(0.2)
        uφ = 5r; p = p0 + T(12.5) * r^2
    elseif r < T(0.4)
        uφ = 2 - 5r; p = p0 + 4 - 4 * log(T(0.2)) + T(12.5) * r^2 - 20r + 4 * log(r)
    else
        uφ = zero(T); p = p0 - 2 + 4 * log(T(2))
    end
    ux = r > 0 ? -uφ * dy / r : zero(T); uy = r > 0 ? uφ * dx / r : zero(T)
    return State{4, T}(one(T), ux, uy, p / (γ - 1) + (ux^2 + uy^2) / 2)
end

struct GreshoIC{T} <: HyperCloud.InitialCondition
    γ::T
    p0::T
end
(ic::GreshoIC{T})(pos) where {T} = gresho_state(T(pos[1]), T(pos[2]), ic.γ, ic.p0)
HyperCloud.build_ic(::Val{:gresho}, conf::Dict, ctx::Dict) = GreshoIC{ctx[:T]}(ctx[:T](conf[:gamma]), ctx[:T](conf[:p0]))
HyperCloud.analytic_closure(eq::HyperCloud.EulerEquation{2, 4, T}, ic::GreshoIC, geom::GeometricDomain) where {T} =
    st -> gresho_state(T(st[1]), T(st[2]), T(eq.gamma), T(ic.p0))

function gresho_config(p0, tmax, Ns, file)
    shared = bench_shared(:tmax => tmax, :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => Ns,
                          :Grid_periodic => true, :IC_name => :gresho, :IC_gamma => 1.4, :IC_p0 => p0)
    return bench_config(shared; reference = :analytical_solution, sources = [joinpath(@__DIR__, "gresho_common.jl"), file])
end
