# Conservation study harness: periodic Burgers / advection on irregular point clouds, built with the
# library builders (so the library scheme is what is measured), plus extra weight functions and helpers
# to extract the MLS divergence coefficients c_ij of the scheme.
#
#   include("study/conservation/Harness.jl"); using .Harness
module Harness

using HyperCloud, StaticArrays, LinearAlgebra, Printf, Random
const HCC = HyperCloud.HyperCloudCore
export setup, mass, evolve!, coefficients, residual_split, WEIGHTS, INTERCEPT, SELF, REG
export Pairs, face_coefficients, rhs!, solve, burgers_exact, burgers_fv_1d, library_rhs
export DISSIPATION, PairsHO, gfd_coefficients, face_coefficients_p, rhs_ho!, solve_ho

# -------------------------------------------------------------------------------------------------
# Extra weight functions (called with the SQUARED distance, like the library ones)
# -------------------------------------------------------------------------------------------------
"exp(-α q²) shifted and rescaled to vanish at the cutoff q = 1 (continuous, kink at q = 1)."
struct ShiftedExp{T} <: HCC.MLSWeightFunction
    alpha::T; range::T; inv_range_sq::T; e::T
end
ShiftedExp(alpha::T, range::T) where {T} = ShiftedExp{T}(alpha, range, inv(range^2), exp(-alpha))
@inline (w::ShiftedExp{T})(d2::Real) where {T} = (q2 = T(d2) * w.inv_range_sq; q2 >= 1 ? zero(T) : (exp(-w.alpha * q2) - w.e) / (1 - w.e))
HCC.get_cutoff(w::ShiftedExp) = w.range

"Wendland C2 kernel (1-q)^4 (4q+1): C2 at the cutoff."
struct WendlandC2{T} <: HCC.MLSWeightFunction
    range::T; inv_range_sq::T
end
WendlandC2(range::T) where {T} = WendlandC2{T}(range, inv(range^2))
@inline (w::WendlandC2{T})(d2::Real) where {T} = (q = sqrt(T(d2) * w.inv_range_sq); q >= 1 ? zero(T) : (1 - q)^4 * (4q + 1))
HCC.get_cutoff(w::WendlandC2) = w.range

"True Gaussian exp(-α q²) truncated at q = 1 (the library :exponential uses a rational approximation)."
struct TruncGauss{T} <: HCC.MLSWeightFunction
    alpha::T; range::T; inv_range_sq::T
end
TruncGauss(alpha::T, range::T) where {T} = TruncGauss{T}(alpha, range, inv(range^2))
@inline (w::TruncGauss{T})(d2::Real) where {T} = (q2 = T(d2) * w.inv_range_sq; q2 > 1 ? zero(T) : exp(-w.alpha * q2))
HCC.get_cutoff(w::TruncGauss) = w.range

_range(c, ctx) = (T = HCC.real_type(ctx); ctx[:interp_range] = T(c[:range]) * ctx[:max_dx]; ctx[:interp_range])
HCC.build_weights(::Val{:shifted_exp}, c::Dict, ctx::Dict) = (T = HCC.real_type(ctx); ShiftedExp(T(c[:alpha]), _range(c, ctx)))
HCC.build_weights(::Val{:wendland}, c::Dict, ctx::Dict) = WendlandC2(_range(c, ctx))
HCC.build_weights(::Val{:gauss}, c::Dict, ctx::Dict) = (T = HCC.real_type(ctx); TruncGauss(T(c[:alpha]), _range(c, ctx)))

const WEIGHTS = (:exponential, :shifted_exp, :wendland, :gauss)

# -------------------------------------------------------------------------------------------------
# Setup of one periodic run on [0,1]^D
# -------------------------------------------------------------------------------------------------
"""
    setup(; D, N, rf, R, weight, alpha, method, mood, mls, pde, ic, volumes, seed, T) -> (eq, pg, ts, params)

`method` is `(stepper, scheme, order)`; `pde` `:burgers` or `:advection` (velocity 1 / (1, 0.5));
`ic(x)` the scalar initial state (default 1 + 0.5 sin(2π x₁) [· sin(2π x₂)]).
"""
function setup(; D = 1, N = 100, rf = 0.2, R = 2.5, weight = :exponential, alpha = 1.0,
               method = (:RK2, :MUSCL, 2), mood = :none, mls = 0, pde = :burgers, ic = nothing,
               volumes = :shepard, seed = 42, T = Float64, limiter = :none)
    p = Dict{Symbol, Any}(last(create_method(method...; limiter = limiter, mood_criterion = mood,
                                             mood_strategy = :EPD1, MOOD_delta_relax = 0.0)))
    merge!(p, Dict{Symbol, Any}(
        :PDE_name => pde === :advection ? :linear : :burgers, :PDE_D => D,
        :PDE_velocities => D == 1 ? ((1.0,),) : ((1.0,), (0.5,)),
        :Grid_domain => :rectangular, :Grid_mins => ntuple(_ -> 0.0, D), :Grid_maxs => ntuple(_ -> 1.0, D),
        :Grid_Ns => ntuple(_ -> N, D), :Grid_periodic => true, :Grid_randomness_factor => ntuple(_ -> rf, D),
        :Grid_seed => seed, :Grid_mover => :none, :Grid_volumes => volumes,
        :Weight_name => weight, :Weight_alpha => alpha, :Weight_range => R,
        :Scheme_MLS_order => mls, :Scheme_upwind_alg_nd => :Classic))
    ctx = Dict{Symbol, Any}(:T => T, :Storage => Array, :ExplicitSources => (), :ImplicitSources => ())
    ns(k) = HyperCloud.extract_namespace(p, k)
    eq = HyperCloud.build_equation(ns(:PDE), ctx); ctx[:Equation] = eq; ctx[:D] = D; ctx[:M] = 1
    ctx[:Domain] = build_domain(ns(:Grid), ctx); ctx[:WeightConf] = ns(:Weight)
    pg = build_particle_grid(ns(:Grid), ctx)
    ctx[:Flux] = build_flux(ns(:Flux), ctx); ctx[:Limiter] = build_limiter(ns(:Limiter), ctx)
    ctx[:MOOD] = build_mood(ns(:MOOD), ctx); ctx[:Scheme] = build_scheme(ns(:Scheme), ctx)
    ctx[:Tableau] = build_tableau(ns(:Time), ctx)
    ts = build_timestepper(ctx; time_conf = ns(:Time))
    f = ic === nothing ? (D == 1 ? x -> 1 + 0.5 * sinpi(2x[1]) : x -> 1 + 0.5 * sinpi(2x[1]) * sinpi(2x[2])) : ic
    HyperCloud.set_states!(x -> State{1, T}(T(f(x))), pg)
    return eq, pg, ts, p
end

"Σ_i V_i U_i (volumes of the grid, or `V` if given) and Σ_i U_i."
mass(pg; V = nothing) = (n = pg.meta.N; vv = V === nothing ? view(pg.core.volumes, 1:n) : V;
                         (sum(vv[i] * pg.rhos[i][1] for i in 1:n), sum(pg.rhos[i][1] for i in 1:n)))

"""
    evolve!(eq, pg, ts; tmax, cfl, V) -> NamedTuple

Runs to `tmax`; returns relative changes of Σ V U (and of Σ U, and of Σ V* U for extra weights `V`),
the number of steps and min/max of the final state.
"""
function evolve!(eq, pg, ts; tmax = 1.0, cfl = 0.2, V = nothing)
    m0, s0 = mass(pg); x0 = V === nothing ? 0.0 : mass(pg; V = V)[1]
    t = 0.0; k = 0
    while t < tmax
        dt = min(cfl * get_time_step(pg, eq, ts.divergence_interpolator), tmax - t)
        dt <= 1e-12 && break
        ts(eq, pg, t, dt); t += dt; k += 1
    end
    m1, s1 = mass(pg); x1 = V === nothing ? 0.0 : mass(pg; V = V)[1]
    u = [pg.rhos[i][1] for i in 1:pg.meta.N]
    return (dVU = (m1 - m0) / m0, dU = (s1 - s0) / s0, dX = V === nothing ? NaN : (x1 - x0) / x0,
            steps = k, t = t, umin = minimum(u), umax = maximum(u))
end

# -------------------------------------------------------------------------------------------------
# MLS divergence coefficients c_ij of the scheme: (∇·F)_i ≈ Σ_j c_ij · (F_j - F_i)
# -------------------------------------------------------------------------------------------------
"""
    coefficients(pg, degree) -> (I, J, C)

Neighbour pairs and the first-derivative MLS coefficients (D-vectors) of the fit of `degree` (the
divergence fit of the scheme; regularised exactly as in the library).
"""
function coefficients(pg::ParticleGrid{D, M, T}, degree::Int; intercept = INTERCEPT[]) where {D, M, T}
    nb = pg.neighbor; n = pg.meta.N
    invL = one(T) / minimum(pg.meta.dx)
    I = Int[]; J = Int[]; C = SVector{D, T}[]
    w0 = nb.weight_func(zero(T))
    for i in 1:n
        sl = nb.ranges[i]
        if intercept
            # least squares WITH a constant term over the neighbours and the particle itself:
            # the derivative weights then sum to zero exactly (Σ_j c_ij = 0)
            Bq = length(HCC.build_basis(Val(degree), zero(SVector{D, T})))
            Ns = zeros(T, Bq + 1, Bq + 1)
            SELF[] && (Ns[1, 1] += w0)
            for k in sl
                p = vcat(one(T), collect(HCC.build_basis(Val(degree), nb.distances[k] * invL))); Ns .+= nb.weights[k] .* (p * p')
            end
            Ns[diagind(Ns)] .+= T(1e-10) * tr(Ns) / (Bq + 1)
            Minv = inv(Symmetric(Ns))
            for k in sl
                p = vcat(one(T), collect(HCC.build_basis(Val(degree), nb.distances[k] * invL)))
                c = Minv * p * nb.weights[k]
                push!(I, i); push!(J, nb.indices[k]); push!(C, SVector{D, T}(ntuple(d -> c[d + 1] * invL, Val(D))))
            end
        else
            Minv, ok = _inv_moment(Val(degree), sl, nb.distances, nb.weights, invL, Val(D), T)
            for k in sl
                p = HCC.build_basis(Val(degree), nb.distances[k] * invL)
                c = Minv * p * nb.weights[k]
                push!(I, i); push!(J, nb.indices[k]); push!(C, SVector{D, T}(ntuple(d -> c[d] * invL, Val(D))))
            end
        end
    end
    return I, J, C
end

"Default of `coefficients(...; intercept)` (the library fit has no constant term)."
const INTERCEPT = Ref(false)
"Relative Tikhonov regularisation of the moment matrices (library: 1e-10)."
const REG = Ref(1e-10)
"With the intercept fit: include the particle itself (then Σ_j c_ij = -c_ii ≠ 0)."
const SELF = Ref(true)
function _inv_moment(::Val{IO}, sl, dists, wts, invL, ::Val{D}, ::Type{T}) where {IO, D, T}
    B = length(HCC.build_basis(Val(IO), zero(SVector{D, T})))
    Ns = zeros(T, B, B)
    for k in sl
        p = HCC.build_basis(Val(IO), dists[k] * invL); Ns .+= wts[k] .* (p * p')
    end
    Ns[diagind(Ns)] .+= T(REG[]) * tr(Ns) / B
    return inv(Symmetric(Ns)), true
end

"""
    residual_split(pg, I, J, C, V) -> NamedTuple

Size of the two conservation defects of `dU_i/dt = -2 Σ_j c_ij·(F*_ij - F_i)` with weights `V`:
`asym` = ‖V_i c_ij + V_j c_ji‖ relative to ‖V_i c_ij‖ (pair antisymmetry) and
`rowsum` = ‖Σ_j c_ij‖ relative to Σ_j ‖c_ij‖ (zero row sums).
"""
function residual_split(pg, I, J, C, V)
    D = length(first(C))
    idx = Dict{Tuple{Int, Int}, Int}((I[k], J[k]) => k for k in eachindex(I))
    num = 0.0; den = 0.0
    for k in eachindex(I)
        kk = get(idx, (J[k], I[k]), 0)
        kk == 0 && continue
        num += norm(V[I[k]] * C[k] + V[J[k]] * C[kk]); den += norm(V[I[k]] * C[k])
    end
    n = pg.meta.N
    rs = [zeros(D) for _ in 1:n]; ra = zeros(n)
    for k in eachindex(I); rs[I[k]] .+= C[k]; ra[I[k]] += norm(C[k]); end
    return (asym = num / den, rowsum = sum(norm.(rs)) / sum(ra))
end

"RHS of the library scheme for the state u (all particles at the maximum order)."
function library_rhs(ts, eq, pg, u)
    n = pg.meta.N; T = eltype(u)
    HCC.update_size!(ts, n, length(pg.neighbor.indices), false); HCC.update_size!(ts.divergence_interpolator, n)
    fill!(view(pg.core.particle_orders, 1:n), HCC._extract_order(ts.divergence_interpolator))
    Y = State{1, T}.(u); K = ts.K_E_stages[1]
    HCC.compute_explicit_rhs!(ts, eq, pg, K, ts.K_X_stages[1], Y, 0.0, true)
    return [K[i][1] for i in 1:n]
end

include("Solver.jl")

end # module
