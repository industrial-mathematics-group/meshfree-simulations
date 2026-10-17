# Forward Euler, RK4MUSCL4(EPD1) + PAD: at the first step that produces a negative pressure, compare the RHS of
# the offending (order-1) particle with that of a pure first-order MUSCL1 scheme on the same state.
using HyperCloud, StaticArrays, Printf, LinearAlgebra
const HCC = HyperCloud.HyperCloudCore
include(joinpath(@__DIR__, "..", "..", "experiments", "ClainShock", "clainshock4.jl"))
cfg = get_config()
function setup(method; T = Float64)
    params = assemble_params(cfg.shared, cfg.methods, Symbol(method))
    params[:Grid_Ns] = (100, 100); params[:Time_stepper] = :Euler
    ctx = Dict{Symbol, Any}(:T => T, :Storage => Array, :ExplicitSources => (), :ImplicitSources => ())
    ns(k) = HyperCloud.extract_namespace(params, k)
    eq = HyperCloud.build_equation(ns(:PDE), ctx); ctx[:Equation] = eq; ctx[:D] = 2; ctx[:M] = 4
    ctx[:Domain] = build_domain(ns(:Grid), ctx); ctx[:WeightConf] = ns(:Weight)
    pg = build_particle_grid(ns(:Grid), ctx)
    ctx[:Flux] = build_flux(ns(:Flux), ctx); ctx[:Limiter] = build_limiter(ns(:Limiter), ctx)
    ctx[:MOOD] = build_mood(ns(:MOOD), ctx); ctx[:Scheme] = build_scheme(ns(:Scheme), ctx)
    ctx[:Tableau] = build_tableau(ns(:Time), ctx)
    ts = build_timestepper(ctx; time_conf = ns(:Time))
    set_initial_conditions!(pg, eq, HyperCloud.build_ic(ns(:IC), ctx))
    return eq, pg, ts
end
eq, pg, ts = setup("RK4MUSCL4(EPD1)")
_, pg1, ts1 = setup("RK2MUSCL1")
pres(s) = (eq.gamma - 1) * (s[4] - (s[2]^2 + s[3]^2) / (2s[1]))
function rhs(ts, pg, Y)
    N = pg.meta.N
    HCC.update_size!(ts, N, length(pg.neighbor.indices), false); HCC.update_size!(ts.divergence_interpolator, N)
    K = ts.K_E_stages[1]
    HCC.compute_explicit_rhs!(ts, eq, pg, K, ts.K_X_stages[1], Y, 0.0, true)
    return copy(view(parent(K), 1:N))
end
function main()
    t = 0.0
    for k in 1:200
        Un = copy(view(pg.rhos, 1:pg.meta.N))
        dt = 0.2 * get_time_step(pg, eq, ts.divergence_interpolator)
        ts(eq, pg, t, dt); t += dt
        bad = [i for i in 1:pg.meta.N if !pg.core.is_boundary[i] && pres(pg.rhos[i]) <= 0]
        isempty(bad) && continue
        i = bad[1]
        orders = copy(pg.core.particle_orders[1:pg.meta.N])
        @printf("step %d: i=%d order %d, p_n=%.4e -> p=%.4e, dt=%.3e\n", k, i, orders[i], pres(Un[i]), pres(pg.rhos[i]), dt)
        # (a) MUSCL4 operator with the final order map, on U_n
        Ka = rhs(ts, pg, Un)
        # (b) pure first order, same state (grid pg1 has the same particles)
        copyto!(pg1.rhos, 1, Un, 1, length(Un)); fill!(pg1.core.particle_orders, 1)
        Kb = rhs(ts1, pg1, Un)
        o_save = copy(pg.core.particle_orders)
        fill!(pg.core.particle_orders, 1)
        Kc = rhs(ts, pg, Un)
        nb = pg.neighbor.indices[pg.neighbor.ranges[i]]
        copyto!(pg.core.particle_orders, o_save)
        for j in nb; pg.core.particle_orders[j] = 1; end           # only i's neighbourhood at order 1
        Kd = rhs(ts, pg, Un)
        copyto!(pg.core.particle_orders, o_save)
        @printf("  K MUSCL4 all order 1: %s\n  K MUSCL4 nbhood order 1: %s\n", Kc[i], Kd[i])
        @printf("  K MUSCL4@order-map: %s\n  K MUSCL1          : %s\n", Ka[i], Kb[i])
        for (m, o) in (("RK2MUSCL2", 2), ("RK4MUSCL3", 3), ("RK4MUSCL4", 4), ("RK4MUSCL5", 5))
            _, pgm, tsm = setup(m)
            copyto!(pgm.rhos, 1, Un, 1, length(Un))
            fill!(pgm.core.particle_orders, 1)
            Km = rhs(tsm, pgm, Un)
            @printf("  K %s all order 1: %s\n", m, Km[i])
        end
        @printf("  FE p with (a): %.4e   with (b): %.4e\n", pres(Un[i] + dt * Ka[i]), pres(Un[i] + dt * Kb[i]))
        @printf("  max dt for (b) positivity check: p(U + 5dt K_b) = %.4e\n", pres(Un[i] + 5dt * Kb[i]))
        return
    end
end
main()
