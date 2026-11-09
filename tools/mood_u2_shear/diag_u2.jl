# Per-test diagnostics of MOOD U2: runs an experiment of the Simulations project (CPU) and, at NSNAP times, forms the
# full-order forward-Euler candidate of the current state and records, for every DMP-violating component of an
# interior particle, the features of the U2 tests. Raw records are serialised for offline evaluation.
#     julia --project=env_headless -t 8 tools/mood_u2_shear/diag_u2.jl CONFIG METHOD TMAX NSNAP OUT.jls [key=value ...]
#     e.g. diag_u2.jl KelvinHelmholtz "RK3MUSCL3(EPD2,AR,MLS2)" 0.3 6 kh.jls "Grid_Ns=(128,128)" T=Float64
using Logging, Printf, Statistics, Serialization, LinearAlgebra
using HyperCloud, StaticArrays
const HCC = HyperCloud.HyperCloudCore
global_logger(ConsoleLogger(stderr, Logging.Warn))

name, method, tend, nsnap, out = ARGS[1], ARGS[2], parse(Float64, ARGS[3]), parse(Int, ARGS[4]), ARGS[5]
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..", "..")
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))
p = HyperCloud.PDEStudioCore.assemble_params(config.shared, config.methods, Symbol(method))
for kv in ARGS[6:end]
    k, v = split(kv, '='; limit = 2)
    p[Symbol(k)] = Core.eval(Main, Meta.parse(v))
end
T = get(p, :T, Float64)
ns(k) = HyperCloud.extract_namespace(p, k)
ctx = Dict{Symbol, Any}(:T => T, :Storage => Array, :ExplicitSources => (), :ImplicitSources => ())
eq = HyperCloud.build_equation(ns(:PDE), ctx); ctx[:Equation] = eq; ctx[:D] = HyperCloud.get_D(eq); ctx[:M] = HyperCloud.get_M(eq)
ctx[:Domain] = HyperCloud.build_domain(ns(:Grid), ctx); ctx[:WeightConf] = ns(:Weight)
pg = HyperCloud.build_run_grid(ns(:Grid), ctx)
ctx[:Flux] = HyperCloud.build_flux(ns(:Flux), ctx); ctx[:Limiter] = HyperCloud.build_limiter(ns(:Limiter), ctx)
ctx[:MOOD] = HyperCloud.build_mood(ns(:MOOD), ctx); ctx[:Scheme] = HyperCloud.build_scheme(ns(:Scheme), ctx)
ctx[:Tableau] = HyperCloud.build_tableau(ns(:Time), ctx)
ts = HyperCloud.build_run_timestepper(ctx; time_conf = ns(:Time))
HyperCloud.set_initial_conditions!(pg, eq, HyperCloud.build_ic(ns(:IC), ctx))
interp = ts.divergence_interpolator
K = p[:Scheme_order]; cfl = T(p[:Time_CFL])
crit = ts.mood.criterion
δ = crit.d
D = HyperCloud.get_D(eq); M = HyperCloud.get_M(eq)

function features(t, dt)
    N = pg.meta.N
    orders = pg.core.particle_orders
    saved = copy(view(orders, 1:N))
    fill!(view(orders, 1:N), K)
    nbs = pg.neighbor
    HCC.update_size!(ts, N, length(nbs.indices), false); HCC.update_size!(interp, N)
    U = copy(view(pg.rhos, 1:N))
    HCC.compute_explicit_rhs!(ts, eq, pg, ts.K_E_stages[1], ts.K_X_stages[1], pg.rhos, t, true)
    KE = parent(ts.K_E_stages[1])
    cand = [U[i] + dt * KE[i] for i in 1:N]
    G = HCC._global_range(HCC.MOODu2(δ; plateau = 1e-3), pg.rhos, N)
    grads = HCC.cached_gradients(interp)
    h = minimum(pg.meta.dx); h2 = h^2
    near_sq = (1.5h)^2
    recs = NamedTuple[]
    nviol_particles = 0; ninterior = 0
    for i in 1:N
        pg.core.is_boundary[i] && continue
        ninterior += 1
        rg = nbs.ranges[i]
        lo, hi = U[i], U[i]
        for k in rg; lo = min.(lo, U[nbs.indices[k]]); hi = max.(hi, U[nbs.indices[k]]); end
        anyv = false
        for m in 1:M
            R = hi[m] - lo[m]
            R >= δ^3 || continue
            v = max(lo[m] - δ - cand[i][m], cand[i][m] - hi[m] - δ)
            v > 0 || continue
            anyv = true
            plateau = R < 1e-3 * G[m]
            # gradient-difference curvature scale over the nearest neighbours (any K ≥ 2)
            gi = SVector{D}(ntuple(d -> grads[i][d][m], D))
            τ1 = zero(T); mis = zero(T); dU = zero(T)
            for k in rg
                x = nbs.distances[k]; r2 = sum(abs2, x)
                j = nbs.indices[k]
                gj = SVector{D}(ntuple(d -> grads[j][d][m], D))
                r2 <= near_sq && (τ1 = max(τ1, h2 * abs(dot(gj - gi, x)) / r2))
                dUj = U[j][m] - U[i][m]
                mis = max(mis, abs(dUj - dot(gi + gj, x) / 2)); dU = max(dU, abs(dUj))
            end
            # cached curvature tests (K ≥ 3)
            τ2 = T(NaN); sign_ok = true; ratio_ok = true; negl = false
            if K >= 3
                sign_ok = true; ratio_ok = true; negl = true; τ2 = zero(T)
                for d in 1:D
                    ci = grads[i][D + d][m]
                    mini = maxi = ci; minA = maxA = abs(ci)
                    for k in rg
                        cj = grads[nbs.indices[k]][D + d][m]
                        mini = min(mini, cj); maxi = max(maxi, cj)
                        if sum(abs2, nbs.distances[k]) <= near_sq
                            minA = min(minA, abs(cj)); maxA = max(maxA, abs(cj))
                        end
                    end
                    τ2 = max(τ2, h2 * max(abs(mini), abs(maxi)))
                    max(abs(mini), abs(maxi)) * h2 < δ && continue
                    negl = false
                    sign_ok &= (mini * h2) * (maxi * h2) > -δ^2
                    ratio_ok &= minA >= crit.ratio * maxA
                end
            end
            push!(recs, (t = Float64(t), i = i, m = m, v = Float64(v), R = Float64(R), G = Float64(G[m]), plateau = plateau,
                         τ1 = Float64(τ1), τ2 = Float64(τ2), sign_ok = sign_ok, ratio_ok = ratio_ok, negl = negl,
                         mis = Float64(mis), dU = Float64(dU), x = Float64.(Tuple(pg.core.positions[i]))))
        end
        anyv && (nviol_particles += 1)
    end
    copyto!(view(orders, 1:N), saved)
    return recs, nviol_particles, ninterior
end

println("$name / $method  N = $(p[:Grid_Ns])  T = $T  K = $K  δ = $δ  up to t = $tend")
allrecs = NamedTuple[]
t = zero(T); k = 0; next = 0.0
while true
    dt = min(cfl * HyperCloud.get_time_step(pg, eq, interp), T(tend) - t)
    if t >= next - 1e-12
        recs, nv, ni = features(t, dt)
        append!(allrecs, recs)
        nr = count(r -> !r.plateau, recs)
        s = K >= 3 ? @sprintf("  sign fail %d  ratio fail %d  both %d  negligible %d", count(r -> !r.plateau && !r.negl && !r.sign_ok, recs),
                              count(r -> !r.plateau && !r.negl && !r.ratio_ok, recs), count(r -> !r.plateau && !r.negl && !r.sign_ok && !r.ratio_ok, recs),
                              count(r -> !r.plateau && r.negl, recs)) : ""
        @printf("t=%.4f step %5d: violating particles %5.2f%%, records %d (non-plateau %d)%s\n", t, k, 100nv / ni, length(recs), nr, s)
        global next += tend / nsnap
    end
    t >= tend - 1e-12 && break
    ts(eq, pg, t, dt); global t += dt; global k += 1
end
serialize(out, (name = name, method = method, K = K, δ = Float64(δ), h = Float64(minimum(pg.meta.dx)), recs = allrecs))
