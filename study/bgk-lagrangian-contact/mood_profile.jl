# MOOD activity along x for a 1-D BGK run: builds the run like HyperCloud's run_bgk_simulation (CPU), steps it by hand
# and records
#   - profiles at t = 0.05, 0.1, 0.2: x, ρ, particle spacing, MOOD particle order o_i and the EPD2 edge order
#     min over N(i) of o (what the edges of i use with EPD2),
#   - the fraction of steps in which a particle is below order K (o_i < K) and in which its EPD2 edge order is below K,
#     accumulated in 200 Eulerian bins of x.
# Output: data/contact/mood_<tag>_profiles.csv and data/contact/mood_<tag>_bins.csv.
#     [SIM_PATCHES=tools/patches/x.jl] julia --project=env_web -t 12 study/bgk-lagrangian-contact/mood_profile.jl CONFIG METHOD TAG [key=value ...]
#     e.g. ... mood_profile.jl bgk_sod_contact EPD2 epd2 Grid_Ns__1=800
using Logging, Printf
using HyperCloud, StaticArrays
global_logger(ConsoleLogger(stderr, Logging.Warn))
const ROOT = abspath(joinpath(@__DIR__, "..", ".."))

for f in split(get(ENV, "SIM_PATCHES", ""), ';'; keepempty = false)   # test patches, as in tools/run_sim.jl
    include(joinpath(ROOT, strip(f)))
end

name, method, tag = ARGS[1], ARGS[2], ARGS[3]
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = ROOT
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))
p = HyperCloud.PDEStudioCore.assemble_params(config.shared, config.methods, Symbol(method))
for (k, v) in config.varied            # first value of every swept key unless overridden below
    p[k] = first(v)
end
for kv in ARGS[4:end]
    k, v = split(kv, '='; limit = 2)
    p[Symbol(k)] = Core.eval(Main, Meta.parse(v))
end
# Grid_Ns__1 style keys set one entry of a tuple parameter
for k in collect(keys(p))
    s = string(k)
    occursin("__", s) || continue
    base, idx = split(s, "__"); base = Symbol(base); idx = parse(Int, idx)
    t = collect(p[base]); t[idx] = p[k]; p[base] = Tuple(t); delete!(p, k)
end
set_backend!(PolyesterBackend())
T = get(p, :T, Float64)
ns(k) = HyperCloud.extract_namespace(p, k)

# --- build (as run_bgk_simulation) ---
ctx = Dict{Symbol, Any}(:T => T, :Storage => Array)
eq = HyperCloud.build_equation(ns(:PDE), ctx)
D = HyperCloud.get_D(eq)
ctx[:Equation] = HyperCloud.euler_limit(eq); ctx[:D] = D; ctx[:M] = D + 2
IC = HyperCloud.build_ic(ns(:IC), ctx)
ctx[:Equation] = eq; ctx[:M] = HyperCloud.get_M(eq)
time_conf = ns(:Time)
ctx[:Tableau] = HyperCloud.build_tableau(time_conf, ctx)
ctx[:ExplicitSources], ctx[:ImplicitSources] = HyperCloud.bgk_sources(eq, ctx[:Tableau])
ctx[:Domain] = HyperCloud.build_domain(ns(:Grid), ctx); ctx[:WeightConf] = ns(:Weight)
pg = HyperCloud.build_run_grid(ns(:Grid), ctx)
ctx[:Flux] = HyperCloud.build_flux(ns(:Flux), ctx); ctx[:Limiter] = HyperCloud.build_limiter(ns(:Limiter), ctx)
ctx[:MOOD] = HyperCloud.build_mood(ns(:MOOD), ctx); ctx[:Scheme] = HyperCloud.build_scheme(ns(:Scheme), ctx)
ts = HyperCloud.build_run_timestepper(ctx; time_conf = time_conf)
HyperCloud.set_initial_conditions!(pg, eq, IC)
div_interp = ts.divergence_interpolator
K = p[:Scheme_order]; cfl = T(p[:Time_CFL]); tmax = T(p[:tmax])
xmin, xmax = p[:Grid_mins][1], p[:Grid_maxs][1]

function nb_min_order(pg, i)
    o = pg.core.particle_orders[i]
    for k in pg.neighbor.ranges[i]
        o = min(o, pg.core.particle_orders[pg.neighbor.indices[k]])
    end
    return o
end

const NB = 200
below_p = zeros(NB); below_e = zeros(NB); visits = zeros(NB)
binof(x) = clamp(1 + floor(Int, (x - xmin) / (xmax - xmin) * NB), 1, NB)
profile_times = (0.05, 0.1, 0.2)
mkpath(joinpath(ROOT, "data", "contact"))
prof = open(joinpath(ROOT, "data", "contact", "mood_$(tag)_profiles.csv"), "w")
println(prof, "t,x,rho,spacing,order,edge_order")

function record_profile(io, pg, t)
    n = pg.meta.N
    inner = [i for i in 1:n if !pg.core.is_boundary[i]]
    xs = [pg.core.positions[i][1] for i in inner]
    perm = sortperm(xs); inner = inner[perm]; xs = xs[perm]
    for (a, i) in enumerate(inner)
        h = a == 1 ? xs[2] - xs[1] : a == length(xs) ? xs[end] - xs[end - 1] : (xs[a + 1] - xs[a - 1]) / 2
        @printf(io, "%.4f,%.6f,%.6f,%.6e,%d,%d\n", t, xs[a], HyperCloud.moments(eq, pg.rhos[i])[1], h,
                pg.core.particle_orders[i], nb_min_order(pg, i))
    end
end

t = zero(T); k = 0; next_prof = 1
println("$name / $method ($tag)  N = $(p[:Grid_Ns])  K = $K  R = $(p[:Weight_range])")
while t < tmax - 1e-12
    dt = min(cfl * HyperCloud.get_time_step(pg, eq, div_interp), tmax - t)
    ts(eq, pg, t, dt); global t += dt; global k += 1
    for i in 1:pg.meta.N
        pg.core.is_boundary[i] && continue
        b = binof(pg.core.positions[i][1])
        visits[b] += 1
        pg.core.particle_orders[i] < K && (below_p[b] += 1)
        nb_min_order(pg, i) < K && (below_e[b] += 1)
    end
    while next_prof <= length(profile_times) && t >= profile_times[next_prof] - 1e-9
        record_profile(prof, pg, profile_times[next_prof]); global next_prof += 1
    end
end
close(prof)
open(joinpath(ROOT, "data", "contact", "mood_$(tag)_bins.csv"), "w") do io
    println(io, "x,below_particle,below_edge,visits")
    for b in 1:NB
        @printf(io, "%.5f,%.5f,%.5f,%d\n", xmin + (b - 0.5) * (xmax - xmin) / NB,
                visits[b] > 0 ? below_p[b] / visits[b] : NaN, visits[b] > 0 ? below_e[b] / visits[b] : NaN, visits[b])
    end
end
inner = count(!, view(pg.core.is_boundary, 1:pg.meta.N))
@printf("steps %d   overall below K: particle %.2f%%  edge %.2f%%\n", k,
        100 * sum(below_p) / sum(visits), 100 * sum(below_e) / sum(visits))
