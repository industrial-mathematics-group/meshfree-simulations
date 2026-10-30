# Sedov blast wave in 2-D (cylindrical; Sedov 1959, Liska & Wendroff 2003): ρ = 1, p = 1e-5 at rest on [-1.2, 1.2]²,
# energy E0 = 0.979264 deposited uniformly (as internal energy) in the disk r < r0 = 0.03, γ = 1.4, t = 1, outflow
# boundaries. Self-similar solution: shock radius R(t) = (E/E0)^(1/4) t^(1/2), i.e. R = 1 at t = 1 for E = E0, density
# jump 6 at the shock, near-vacuum and very high temperature at the centre: a positivity test.
# Post-processing statistic (sedov_stats!) at t_end, stored as one 4-vector :sedov_shock = (R_mean, R_pred, spread,
# ρ_max): R_mean the mean over 64 angular sectors of the radius of the density maximum, R_pred = (E_dep/E0)^(1/4) t^(1/2)
# with the energy E_dep actually deposited on the particles (Σ h² (E_i − E_ambient) at t = 0), spread =
# (max − min)/mean of the sector radii (asymmetry of the shock), ρ_max the largest density.
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

const SEDOV_E0 = 0.979264      # γ = 1.4, ρ = 1: R = 1 at t = 1

struct SedovIC{T} <: HyperCloud.InitialCondition
    r0::T
    E_in::T       # energy density inside r0 (blast + ambient)
    E_out::T      # ambient energy density p / (γ − 1)
end
(ic::SedovIC{T})(pos) where {T} = State{4, T}(one(T), zero(T), zero(T), pos[1]^2 + pos[2]^2 < ic.r0^2 ? ic.E_in : ic.E_out)
function HyperCloud.build_ic(::Val{:sedov}, conf::Dict, ctx::Dict)
    T = ctx[:T]; γ, r0, p = conf[:gamma], conf[:r0], conf[:p_ambient]
    return SedovIC{T}(T(r0), T(conf[:E0] / (π * r0^2) + p / (γ - 1)), T(p / (γ - 1)))
end

function sedov_stats!(data::HyperCloud.AbstractSimData{D, DS, M, T}) where {D, DS, M, T}
    haskey(data.stats, :sedov_shock) && return false
    p = data.params; γ = p[:PDE_gamma]
    h2 = prod(Float64.(data.domain.spacing[1:DS]))
    E_dep = sum(u -> Float64(u[4]) - p[:IC_p_ambient] / (γ - 1), data.u[1]) * h2
    t = Float64(data.t[end])
    R_pred = (E_dep / SEDOV_E0)^(1 / 4) * sqrt(t)
    nb = 64
    rmax = zeros(nb); ρmax = zeros(nb)
    for (x, u) in zip(data.x[end], data.u[end])
        b = clamp(floor(Int, (atan(Float64(x[2]), Float64(x[1])) + π) / (2π) * nb) + 1, 1, nb)
        if u[1] > ρmax[b]
            ρmax[b] = u[1]; rmax[b] = hypot(Float64(x[1]), Float64(x[2]))
        end
    end
    R = sum(rmax) / nb
    add_stat!(data, :sedov_shock, SVector{M, T}(R, R_pred, (maximum(rmax) - minimum(rmax)) / R, maximum(ρmax)), Symbol[])
    return true
end

function get_config()
    shared = bench_shared(:tmax => 1.0, :Grid_mins => (-1.2, -1.2), :Grid_maxs => (1.2, 1.2), :Grid_Ns => (400, 400),
                          :Grid_bc_map => bench_outflow(), :IC_name => :sedov, :IC_gamma => 1.4, :IC_r0 => 0.03,
                          :IC_E0 => SEDOV_E0, :IC_p_ambient => 1e-5)
    return bench_config(shared; postprocess = :sedov_stats!, sources = [@__FILE__])
end
