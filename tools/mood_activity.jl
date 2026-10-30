# MOOD activity of a run: steps an experiment (built like run_direct_simulation, CPU) and reports, at regular intervals,
# the fraction of interior particles whose MOOD order after the step is below the scheme order K and the fraction at
# order 1, with the mean and quartiles of a position function of the dropped particles (to see where MOOD acts).
#     julia --project=env_headless -t 12 tools/mood_activity.jl CONFIG METHOD TMAX [key=value ...]
#     e.g. tools/mood_activity.jl KelvinHelmholtz "RK2MUSCL2(EPD2,AR,MLS2)" 0.5 "Grid_Ns=(256,256)" T=Float64
# env MOOD_WHERE: Julia expression in x, y for the position statistic (default: distance to the KH interfaces).
using Logging, Printf, Statistics
using HyperCloud, StaticArrays
const HCC = HyperCloud.HyperCloudCore
global_logger(ConsoleLogger(stderr, Logging.Warn))

name, method, tend = ARGS[1], ARGS[2], parse(Float64, ARGS[3])
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))
p = HyperCloud.PDEStudioCore.assemble_params(config.shared, config.methods, Symbol(method))
for kv in ARGS[4:end]
    k, v = split(kv, '='; limit = 2)
    p[Symbol(k)] = Core.eval(Main, Meta.parse(v))
end
T = get(p, :T, Float64)
where_f = Core.eval(Main, Meta.parse("(x, y) -> " * get(ENV, "MOOD_WHERE", "min(abs(y - 0.25), abs(y - 0.75))")))

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
div_interp = ts.divergence_interpolator
K = p[:Scheme_order]; cfl = T(p[:Time_CFL])

println("$name / $method  N = $(p[:Grid_Ns])  T = $T  K = $K  up to t = $tend")
println("    t       steps  below K   order 1   position of the dropped particles: mean / q25 / q50 / q75")
t = zero(T); k = 0; next = 0.0
while t < tend
    dt = min(cfl * HyperCloud.get_time_step(pg, eq, div_interp), T(tend) - t)
    ts(eq, pg, t, dt); global t += dt; global k += 1
    if t >= next - 1e-12 || t >= tend
        n = pg.meta.N
        ib = view(pg.core.is_boundary, 1:n); o = view(pg.core.particle_orders, 1:n); x = view(pg.core.positions, 1:n)
        inner = findall(!, ib)
        low = [i for i in inner if o[i] < K]
        w = [where_f(Float64(x[i][1]), Float64(x[i][2])) for i in low]
        q = isempty(w) ? (NaN, NaN, NaN) : (quantile(w, 0.25), quantile(w, 0.5), quantile(w, 0.75))
        @printf("%8.4f  %6d  %7.2f%%  %7.2f%%   %.4f / %.4f / %.4f / %.4f\n", t, k, 100 * length(low) / length(inner),
                100 * count(i -> o[i] == 1, inner) / length(inner), isempty(w) ? NaN : mean(w), q...)
        global next += tend / 10
    end
end
