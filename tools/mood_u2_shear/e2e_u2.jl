# End-to-end check of MOOD U2 variants on Simulations experiments (CPU): metrics at TMAX plus MOOD activity.
#     julia --project=env_headless -t 8 tools/mood_u2_shear/e2e_u2.jl CONFIG METHOD TMAX [key=value ...]
using Logging, Printf, Statistics, LinearAlgebra
using HyperCloud, StaticArrays
const HCC = HyperCloud.HyperCloudCore
global_logger(ConsoleLogger(stderr, Logging.Warn))

name, method, tend = ARGS[1], ARGS[2], parse(Float64, ARGS[3])
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..", "..")
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))
p = HyperCloud.PDEStudioCore.assemble_params(config.shared, config.methods, Symbol(method))
extra = String[]
for kv in ARGS[4:end]
    k, v = split(kv, '='; limit = 2)
    p[Symbol(k)] = Core.eval(Main, Meta.parse(v)); push!(extra, kv)
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
K = p[:Scheme_order]; cfl = T(p[:Time_CFL]); D = HyperCloud.get_D(eq)
N0 = pg.meta.N
U0 = copy(view(pg.rhos, 1:N0)); X0 = copy(view(pg.core.positions, 1:N0))
γ = Float64(get(p, :PDE_gamma, 1.4))
pres(U) = D == 1 ? (γ - 1) * (U[3] - U[2]^2 / (2U[1])) : (γ - 1) * (U[4] - (U[2]^2 + U[3]^2) / (2U[1]))

t = zero(T); k = 0; below = Float64[]; first = Float64[]
t0 = time()
while t < tend - 1e-12
    dt = min(cfl * HyperCloud.get_time_step(pg, eq, interp), T(tend) - t)
    ts(eq, pg, t, dt); global t += dt; global k += 1
    if k % 5 == 0
        n = pg.meta.N; o = view(pg.core.particle_orders, 1:n); ib = view(pg.core.is_boundary, 1:n)
        ni = count(!, ib)
        push!(below, count(i -> !ib[i] && o[i] < K, 1:n) / ni); push!(first, count(i -> !ib[i] && o[i] == 1, 1:n) / ni)
    end
end
wall = time() - t0
n = pg.meta.N; U = view(pg.rhos, 1:n); X = view(pg.core.positions, 1:n)
ρ = [Float64(u[1]) for u in U]
function kh_amp(X, U)
    s = c = d = 0.0
    for (x, u) in zip(X, U)
        Xx, Y = Float64(x[1]), Float64(x[2]); w = exp(-4π * abs(Y - (Y < 0.5 ? 0.25 : 0.75)))
        v = Float64(u[3]) / Float64(u[1]); s += v * sinpi(4Xx) * w; c += v * cospi(4Xx) * w; d += w
    end
    return 2 * sqrt((s / d)^2 + (c / d)^2)
end
metric = if name == "KelvinHelmholtz"
    @sprintf("KH mode %.4e", kh_amp(X, U))
elseif startswith(name, "Gresho")
    @sprintf("L1(ρu) %.4e  L1(ρ) %.4e", mean(abs(Float64(U[i][2] - U0[i][2])) + abs(Float64(U[i][3] - U0[i][3])) for i in 1:n) / 2 * 2,
             mean(abs(Float64(U[i][1] - U0[i][1])) for i in 1:n))
elseif name == "Sod1D"
    ex = analytical_solution(p)
    err = mean(abs(ex(SVector(Tuple(X[i])..., Float64(t)))[1] - ρ[i]) for i in 1:n)
    @sprintf("L1(ρ) %.4e  over %.2e  under %.2e", err, maximum(ρ) - 1, 0.125 - minimum(ρ))
elseif name == "Implosion"
    key(x) = (round(Int, Float64(x[1]) * 1e6), round(Int, Float64(x[2]) * 1e6))
    idx = Dict(key(X[i]) => i for i in 1:n)
    sym = mean(abs(ρ[i] - ρ[get(idx, (key(X[i])[2], key(X[i])[1]), i)]) for i in 1:n)
    @sprintf("sym %.3e", sym)
else
    ex = try analytical_solution(p) catch; nothing end
    ex === nothing ? "" : @sprintf("L1(ρ) %.4e", mean(abs(ex(SVector(Tuple(X[i])..., Float64(t)))[1] - ρ[i]) for i in 1:n))
end
lo0, hi0 = extrema(Float64(u[1]) for u in U0)
imin, imax = argmin(ρ), argmax(ρ)
@printf("  argmin ρ at %s (boundary %s), argmax ρ at %s; particles with ρ < init min: %d\n", string(round.(Float64.(Tuple(X[imin])), digits = 3)),
        pg.core.is_boundary[imin], string(round.(Float64.(Tuple(X[imax])), digits = 3)), count(<(lo0 - 1e-6), ρ))
@printf("%-16s %-26s %-18s t=%.3f steps=%5d: %s  ρ∈[%.4f, %.4f] (init [%.4f, %.4f])  min p %.3e  below K %.2f%%  order 1 %.2f%%  wall %.0fs\n",
        name, method, join(extra, " "), t, k, metric, minimum(ρ), maximum(ρ), lo0, hi0, minimum(pres, U),
        100mean(below), 100mean(first), wall)
