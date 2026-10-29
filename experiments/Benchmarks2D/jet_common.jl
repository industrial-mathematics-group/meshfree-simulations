# High-Mach astrophysical jets (Ha, Gardner et al. 2005; Zhang & Shu 2010, positivity-preserving schemes): γ = 5/3,
# ambient gas (ρ, u, v, p) = (0.5, 0, 0, 0.4127) at rest, a jet (5, u_jet, 0, 0.4127) enters through the part
# |y| < 0.05 of the left boundary, outflow everywhere else. Kinetic energy dwarfs the internal energy (Mach 80 / 2000
# relative to the ambient sound speed), so the pressure p = (γ − 1)(E − ρ|u|²/2) is a small difference of large
# numbers: a positivity test, in particular in single precision. Configs: MachJet80.jl, MachJet2000.jl.
# The jet state also fills a strip x < 0.005 (2 … 3 spacings) of the domain at t = 0: the CFL estimate of
# HyperCloudCore 0.7.8 uses only each interior particle's own state (boundary particles are skipped), so with ambient
# gas at rest everywhere inside, the first step would be ~25× (Mach 80) / ~700× (Mach 2000) too large for the jet
# entering through the boundary in that step.
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

"Left boundary: zero-gradient outflow, except the jet state for |y| < half_width."
struct JetInflowBC <: AbstractBoundaryCondition
    jet::NTuple{4, Float64}
    half_width::Float64
end
@inline function _jet_bc!(i, buf, pos, tags, is_boundary, tag, jet, hw)
    if is_boundary[i] && tags[i] == tag && abs(pos[i][2]) < hw
        buf[i] = jet
    end
    return nothing
end
function (bc::JetInflowBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t)
    S = eltype(rhos_buffer); T = eltype(S)
    foreach_particle(_jet_bc!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.positions, pg.core.tags, pg.core.is_boundary, tag, S(bc.jet), T(bc.half_width)))
    return nothing
end

function jet_config(u_jet, mins, maxs, Ns, tmax, file)
    γ = 5 / 3
    ambient = bench_cons(γ, 0.5, 0, 0, 0.4127)
    jet = bench_cons(γ, 5.0, u_jet, 0, 0.4127)
    bc = Dict{Int, AbstractBoundaryCondition}(1 => JetInflowBC(jet, 0.05),
                                              2 => OutflowBC(), 3 => OutflowBC(), 4 => OutflowBC())
    shared = bench_shared(:PDE_gamma => γ, :tmax => tmax, :Grid_mins => mins, :Grid_maxs => maxs, :Grid_Ns => Ns,
                          :Grid_bc_map => bc, :IC_name => :box, :IC_u_bg => ambient, :IC_u_box => jet,
                          :IC_mins => (mins[1] - 1.0, -0.05), :IC_maxs => (mins[1] + 0.005, 0.05))
    return bench_config(shared; sources = [joinpath(@__DIR__, "jet_common.jl"), file])
end
