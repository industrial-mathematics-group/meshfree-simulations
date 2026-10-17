# Steps clainshock4 (100x100, RK4/RK3 MUSCL4 EPD1, PAD branch) and reports the first particles with negative
# pressure: their MOOD order and the orders of their neighbours.
using HyperCloud, StaticArrays, Printf
const HCC = HyperCloud.HyperCloudCore
include(joinpath(@__DIR__, "..", "..", "experiments", "ClainShock", "clainshock4.jl"))
cfg = get_config()
params = assemble_params(cfg.shared, cfg.methods, Symbol(get(ENV, "METHOD", "RK4MUSCL4(EPD1)")))
params[:Weight_range] = parse(Float64, get(ENV, "RANGE", "5.5"))
params[:Scheme_MLS_order] = parse(Int, get(ENV, "DIVORDER", "4"))
params[:Grid_Ns] = (100, 100); params[:Time_CFL] = 0.2
params[:Time_stepper] = Symbol(get(ENV, "STEPPER", "RK4"))
params[:Time_restart_on_drop] = get(ENV, "RESTART", "0") == "1"
T = Float64
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
γ = eq.gamma
pres(s) = (γ - 1) * (s[4] - (s[2]^2 + s[3]^2) / (2s[1]))
function run!(ts, eq, pg)
t = 0.0; k = 0
while t < 0.3
    dt = min(0.2 * get_time_step(pg, eq, ts.divergence_interpolator), 0.3 - t)
    ts(eq, pg, t, dt); t += dt; k += 1
    N = pg.meta.N
    bad = [i for i in 1:N if !pg.core.is_boundary[i] && pres(pg.rhos[i]) <= 0]
    isempty(bad) && continue
    @printf("step %d t=%.4f: %d particles with p<=0\n", k, t, length(bad))
    for i in bad[1:min(end, 5)]
        nb = pg.neighbor.indices[pg.neighbor.ranges[i]]
        ords = pg.core.particle_orders[nb]
        @printf("  i=%d x=(%.3f,%.3f) p=%.3e order=%d  nb orders: min %d max %d, #order1 %d/%d, boundary nbs %d\n",
                i, pg.core.positions[i]..., pres(pg.rhos[i]), pg.core.particle_orders[i], minimum(ords), maximum(ords),
                count(==(1), ords), length(nb), count(j -> pg.core.is_boundary[j], nb))
    end
    k >= parse(Int, get(ENV, "MAXSTEPS", "1000")) && break
    get(ENV, "STOP_FIRST", "1") == "1" && break
end
println("done at step $k t=$t")
end
run!(ts, eq, pg)
