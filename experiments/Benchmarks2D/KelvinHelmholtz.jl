# Kelvin–Helmholtz instability, the smooth, well-posed setup of McNally, Lyra & Passy (2012): periodic [0, 1]², γ = 5/3,
# p = 2.5, a dense band (ρ = 2, u = −0.5) in 0.25 < y < 0.75 between light gas (ρ = 1, u = 0.5), interfaces smoothed
# with width L = 0.025, seeded with v = 0.01 sin(4πx); t = 2. Numerical dissipation sets how fast the seeded mode grows
# and how fine the roll-ups are (no pointwise grid convergence late in time, so compare statistics).
# Post-processing statistic (kh_mode!): the amplitude of the seeded y-velocity mode (McNally et al. eq. 10–13),
#   :kh_mode  time series  M(t) = 2 sqrt((s/d)² + (c/d)²),  s, c, d = Σ_i h² {v_i sin(4πx_i), v_i cos(4πx_i), 1} e^(−4π|y_i − y_I|),
# with y_I = 0.25 for y < 0.5 and 0.75 above; :final_kh_mode its value at t_end. M(0) = 0.01.
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

@inline function kh_state(x::T, y::T, γ::T) where {T}
    ρ1, ρ2, u1, u2, L = one(T), T(2), T(0.5), T(-0.5), T(0.025)
    ρm, um = (ρ1 - ρ2) / 2, (u1 - u2) / 2
    if y < T(0.25)
        e = exp((y - T(0.25)) / L); ρ = ρ1 - ρm * e; u = u1 - um * e
    elseif y < T(0.5)
        e = exp((T(0.25) - y) / L); ρ = ρ2 + ρm * e; u = u2 + um * e
    elseif y < T(0.75)
        e = exp((y - T(0.75)) / L); ρ = ρ2 + ρm * e; u = u2 + um * e
    else
        e = exp((T(0.75) - y) / L); ρ = ρ1 - ρm * e; u = u1 - um * e
    end
    v = T(0.01) * sin(4 * T(π) * x)
    return State{4, T}(ρ, ρ * u, ρ * v, T(2.5) / (γ - 1) + ρ * (u^2 + v^2) / 2)
end

struct KelvinHelmholtzIC{T} <: HyperCloud.InitialCondition
    γ::T
end
(ic::KelvinHelmholtzIC{T})(pos) where {T} = kh_state(T(pos[1]), T(pos[2]), ic.γ)
HyperCloud.build_ic(::Val{:kelvin_helmholtz}, conf::Dict, ctx::Dict) = KelvinHelmholtzIC{ctx[:T]}(ctx[:T](conf[:gamma]))

function kh_mode!(data::HyperCloud.AbstractSimData{D, DS, M, T}) where {D, DS, M, T}
    haskey(data.stats, :final_kh_mode) && return false
    amp = map(eachindex(data.t)) do n
        s = c = d = 0.0
        for (x, U) in zip(data.x[n], data.u[n])
            X, Y = Float64(x[1]), Float64(x[2])
            w = exp(-4π * abs(Y - (Y < 0.5 ? 0.25 : 0.75)))
            v = Float64(U[3]) / Float64(U[1])
            s += v * sinpi(4X) * w; c += v * cospi(4X) * w; d += w
        end
        T(2 * sqrt((s / d)^2 + (c / d)^2))
    end
    add_stat!(data, :kh_mode, [SVector{M, T}(ntuple(_ -> a, M)) for a in amp], :time)
    add_stat!(data, :final_kh_mode, SVector{M, T}(ntuple(_ -> last(amp), M)), Symbol[])
    return true
end

function get_config()
    shared = bench_shared(:PDE_gamma => 5 / 3, :tmax => 2.0, :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0),
                          :Grid_Ns => (512, 512), :Grid_periodic => true, :IC_name => :kelvin_helmholtz, :IC_gamma => 5 / 3)
    return bench_config(shared; postprocess = :kh_mode!, sources = [@__FILE__])
end
