# Discrete conservation Σ V_i U_i of a periodic, jittered 2-D advection of a box (MOOD active) for a fixed
# MLS divergence order vs the coupled one (Scheme_MLS_order = 0).
#   julia --project=env_pad tools/debug/conservation.jl
using HyperCloud, StaticArrays, Printf
const HCC = HyperCloud.HyperCloudCore
function run(method, mls; N = 100, R = 5.5, rf = 0.2, tmax = 0.5, T = Float64, mood = :U2, alpha = 1.0)
    p = Dict{Symbol, Any}(last(create_method(method...; mood_criterion = mood, mood_strategy = :EPD1, MOOD_delta_relax = 0.0)))
    merge!(p, Dict{Symbol, Any}(:PDE_name => :linear, :PDE_D => 2, :PDE_velocities => ((1.0,), (0.5,)),
        :Grid_domain => :rectangular, :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => (N, N),
        :Grid_periodic => true, :Grid_randomness_factor => (rf, rf), :Grid_seed => 42, :Grid_mover => :none,
        :Weight_name => :exponential, :Weight_alpha => alpha, :Weight_range => R,
        :IC_name => :box, :IC_u_bg => (0.0,), :IC_u_box => (1.0,), :IC_mins => (0.25, 0.25), :IC_maxs => (0.6, 0.6),
        :Scheme_MLS_order => mls, :Scheme_upwind_alg_nd => :Classic))
    ctx = Dict{Symbol, Any}(:T => T, :Storage => Array, :ExplicitSources => (), :ImplicitSources => ())
    ns(k) = HyperCloud.extract_namespace(p, k)
    eq = HyperCloud.build_equation(ns(:PDE), ctx); ctx[:Equation] = eq; ctx[:D] = 2; ctx[:M] = 1
    ctx[:Domain] = build_domain(ns(:Grid), ctx); ctx[:WeightConf] = ns(:Weight)
    pg = build_particle_grid(ns(:Grid), ctx)
    ctx[:Flux] = build_flux(ns(:Flux), ctx); ctx[:Limiter] = build_limiter(ns(:Limiter), ctx)
    ctx[:MOOD] = build_mood(ns(:MOOD), ctx); ctx[:Scheme] = build_scheme(ns(:Scheme), ctx)
    ctx[:Tableau] = build_tableau(ns(:Time), ctx)
    ts = build_timestepper(ctx; time_conf = ns(:Time))
    set_initial_conditions!(pg, eq, HyperCloud.build_ic(ns(:IC), ctx))
    n = pg.meta.N
    if get(ENV, "RESID", "0") == "1"
        HCC.update_size!(ts, n, length(pg.neighbor.indices), false); HCC.update_size!(ts.divergence_interpolator, n)
        fill!(pg.core.particle_orders, 1)
        for (name, f) in (("sin", x -> sin(2pi * x[1]) * cos(4pi * x[2])), ("box", x -> pg.rhos[1][1] * 0 + (0.25 <= x[1] <= 0.6 && 0.25 <= x[2] <= 0.6)))
            Y = [SVector(Float64(f(pg.core.positions[i]))) for i in 1:n]
            K = ts.K_E_stages[1]
            HCC.compute_explicit_rhs!(ts, eq, pg, K, ts.K_X_stages[1], Y, 0.0, true)
            Kv = [K[i][1] for i in 1:n]
            println("N=$n rf=$rf alpha=$alpha $name: Σ V_i K_i = $(sum(pg.core.volumes[1:n] .* Kv)),  Σ|V_i K_i| = $(sum(abs.(pg.core.volumes[1:n] .* Kv)))")
        end
        return
    end
    if get(ENV, "GEOM", "0") == "1"
        nb = pg.neighbor; pos = pg.core.positions[1:n]
        sets = [Set(nb.indices[nb.ranges[i]]) for i in 1:n]
        asym = count(i -> any(j -> !(i in sets[j]), sets[i]), 1:n)
        outside = count(x -> any(x .< 0) || any(x .>= 1), pos)
        dmis = 0
        for i in 1:n, k in nb.ranges[i]
            j = nb.indices[k]; kk = findfirst(==(i), nb.indices[nb.ranges[j]])
            kk === nothing && continue
            dji = nb.distances[nb.ranges[j][kk]]
            maximum(abs.(nb.distances[k] .+ dji)) > 1e-12 && (dmis += 1)
        end
        println("N=$n rf=$rf: particles with asymmetric neighbours $asym, outside [0,1) $outside, pairs with d_ij != -d_ji $dmis, i in own list $(count(i -> i in sets[i], 1:n))")
        return
    end
    mass() = sum(pg.core.volumes[i] * pg.rhos[i][1] for i in 1:n)
    plain() = sum(pg.rhos[i][1] for i in 1:n)
    m0 = mass(); s0 = plain(); v0 = copy(pg.core.volumes[1:n]); t = 0.0; k = 0; dropped = 0
    while t < tmax
        dt = min(0.2 * get_time_step(pg, eq, ts.divergence_interpolator), tmax - t)
        ts(eq, pg, t, dt); t += dt; k += 1
        dropped += count(HCC.mood_flags(ts))
    end
    @printf("N=%d rf=%.2f mood=%-4s %-8s MLS_order=%d  steps=%4d  MOOD drops/step=%7.1f  rel. ΣVU change=%+.3e  ΣU change=%+.3e  V changed=%s  min=%+.3e max=%.4f\n",
            N, rf, mood, join(string.(method[2:3])), mls, k, dropped / k, (mass() - m0) / m0, (plain() - s0) / s0, pg.core.volumes[1:n] != v0,
            minimum(pg.rhos[i][1] for i in 1:n), maximum(pg.rhos[i][1] for i in 1:n))
end
for alpha in (1.0, 4.0, 10.0), rf in (0.01, 0.2)
    run((:RK2, :MUSCL, 1), 0; rf = rf, mood = :none, N = 50, alpha = alpha)
end
