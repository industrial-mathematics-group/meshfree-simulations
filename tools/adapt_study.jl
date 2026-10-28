# Particle-management study: runs one experiment/method with the run pipeline's builders, but with its own time
# loop, so that every particle insertion/merge is recorded with the remap's own mass residual. Per step it logs
#   t, N, inserted, removed, Δm_remap (adapt_particles!'s mass_residual = change of Σ V_i U_i over the call),
#   Δm_step (change of Σ V_i U_i over the time step), Σ V_i (total volume)
# so that the conservation defect of the remap and of the scheme are separated. Writes
# <out>/<tag>.csv (per step) and <out>/<tag>.jld2 (summary + snapshots), prints RESULT lines.
#
#   julia --project=env_headless tools/adapt_study.jl NAME METHOD TAG [key=value ...]
#     NAME   experiment (as for @sim), METHOD a method label of its config, TAG output name
#     key=value overrides of the merged parameters (Julia literals), e.g. Grid_Ns=(400,) Grid_mover=:lagrangian
#       Adapt_enabled=true Adapt_max_gap=1.6 Adapt_min_gap=0.5 Adapt_samples_per_dx=2 Adapt_exclude_reduced_order=true
#     env ADAPT_OUT (default results/adapt), ADAPT_SNAPSHOTS (default 10)
using HyperCloud, StaticArrays, Printf, Logging
const HCC = HyperCloud.HyperCloudCore
const JLD2 = HyperCloud.PDEStudioCore.JLD2
global_logger(ConsoleLogger(stderr, Logging.Warn))
get(ENV, "ADAPT_BGK_ALE_FLUX", "0") == "1" && include(joinpath(@__DIR__, "patches", "bgk_ale_flux.jl"))

name, method, tag = ARGS[1], ARGS[2], ARGS[3]
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa Function && (config = Base.invokelatest(config))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))
p = HyperCloud.PDEStudioCore.assemble_params(config.shared, config.methods, Symbol(method))
for kv in ARGS[4:end]
    k, v = split(kv, '='; limit = 2)
    p[Symbol(k)] = Core.eval(Main, Meta.parse(v))
end
T = get(p, :T, Float64)
out = get(ENV, "ADAPT_OUT", joinpath(@__DIR__, "..", "results", "adapt")); mkpath(out)

ns(k) = HyperCloud.extract_namespace(p, k)
ctx = Dict{Symbol, Any}(:T => T, :Storage => Array, :ExplicitSources => (), :ImplicitSources => ())
eq = HyperCloud.build_equation(ns(:PDE), ctx); ctx[:Equation] = eq; ctx[:D] = HyperCloud.get_D(eq); ctx[:M] = HyperCloud.get_M(eq)
const IS_BGK = p[:PDE_name] == :bgk
if IS_BGK   # as run_bgk_simulation: the IC is macroscopic (built with the Euler limit), the collision a source term
    ctx[:Equation] = HyperCloud.euler_limit(eq); ctx[:M] = ctx[:D] + 2
    IC0 = HyperCloud.build_ic(ns(:IC), ctx)
    ctx[:Equation] = eq; ctx[:M] = HyperCloud.get_M(eq)
    ctx[:Tableau] = HyperCloud.build_tableau(ns(:Time), ctx)
    ctx[:ExplicitSources], ctx[:ImplicitSources] = HyperCloud.bgk_sources(eq, ctx[:Tableau])
end
ctx[:Domain] = HyperCloud.build_domain(ns(:Grid), ctx); ctx[:WeightConf] = ns(:Weight)
pg = HyperCloud.build_run_grid(ns(:Grid), ctx)
ad = HyperCloud.build_adaptivity(ns(:Adapt), ctx)
ctx[:Flux] = HyperCloud.build_flux(ns(:Flux), ctx); ctx[:Limiter] = HyperCloud.build_limiter(ns(:Limiter), ctx)
ctx[:MOOD] = HyperCloud.build_mood(ns(:MOOD), ctx); ctx[:Scheme] = HyperCloud.build_scheme(ns(:Scheme), ctx)
IS_BGK || (ctx[:Tableau] = HyperCloud.build_tableau(ns(:Time), ctx))
# The run pipeline does not pass the grid's volume estimator to the stepper (nothing sets context[:VolumeEstimator],
# so `Grid_volumes` only affects the initial volumes, HyperCloudCore 0.7.3); do it here so the option applies.
pkgversion(HCC) < v"0.7.5" && (ctx[:VolumeEstimator] = HyperCloud.build_volume_estimator(ns(:Grid), ctx))   # fixed in 0.7.5
ts = HyperCloud.build_run_timestepper(ctx; time_conf = ns(:Time))
HyperCloud.set_initial_conditions!(pg, eq, IS_BGK ? IC0 : HyperCloud.build_ic(ns(:IC), ctx))
div_interp = ts.divergence_interpolator
order = HCC._extract_order(div_interp)

tmax = T(p[:tmax]); cfl = T(p[:Time_CFL]); nsnap = parse(Int, get(ENV, "ADAPT_SNAPSHOTS", "10"))
interior(pg) = .!collect(view(pg.core.is_boundary, 1:pg.meta.N))
# conserved content; for BGK the moments (mass, momentum, energy) of Σ V_i f_i (moments are linear in f)
mass(pg) = IS_BGK ? HyperCloud.moments(eq, HyperCloud.total_mass(pg)) : HyperCloud.total_mass(pg)
vol(pg) = sum(view(pg.core.volumes, 1:pg.meta.N))
m0 = mass(pg); N0 = pg.meta.N
mint0 = let n = pg.meta.N, ib = collect(view(pg.core.is_boundary, 1:n)), V = collect(view(pg.core.volumes, 1:n)), F = collect(view(pg.rhos, 1:n))
    mi = sum(V[i] * F[i] for i in 1:n if !ib[i]); IS_BGK ? HyperCloud.moments(eq, mi) : mi
end
io = open(joinpath(out, tag * ".csv"), "w")
println(io, "step,t,N,inserted,removed,dm_remap,dm_step,volume,mass")
snapstate(pg) = IS_BGK ? [HyperCloud.moments(eq, u) for u in view(pg.rhos, 1:pg.meta.N)] : collect(view(pg.rhos, 1:pg.meta.N))
snaps = Any[(0.0, collect(view(pg.core.positions, 1:N0)), snapstate(pg))]
sum_revol = zero(m0); t = zero(T); k = 0; ins = 0; rem = 0; sum_remap = zero(m0); sum_step = zero(m0); max_remap = zero(T)
next_snap = tmax / nsnap
status = "ok"
elapsed = @elapsed while t < tmax
    dt = min(cfl * HyperCloud.get_time_step(pg, eq, div_interp), tmax - t)
    dt <= T(1e-12) && (global status = "dt collapse"; break)
    mb = mass(pg)
    ts(eq, pg, t, dt)
    ma = mass(pg)
    evs = (ad !== nothing && get(ENV, "ADAPT_EVENTLOG", "0") == "1") ? HCC._detect_events(pg, ad, order) : nothing
    pg_before = evs !== nothing && haskey(ENV, "ADAPT_HOLECHECK") ? deepcopy(pg) : nothing
    r = ad === nothing ? (inserted = 0, removed = 0, mass_residual = zero(m0)) : HyperCloud.adapt_particles!(pg, ad, order)
    IS_BGK && ad !== nothing && (r = (inserted = r.inserted, removed = r.removed, mass_residual = HyperCloud.moments(eq, r.mass_residual)))
    if evs !== nothing && maximum(abs, r.mass_residual) > 1e-11 * maximum(abs, m0) && haskey(ENV, "ADAPT_HOLECHECK")
        # coverage around each event of the bad call, on the cloud BEFORE the call (positions saved below)
        Rw = p[:Weight_range] * minimum(pg.meta.dx); idx = Vector{Int}(undef, 4096); dv = Vector{Space{ctx[:D], T}}(undef, 4096)
        for ev in evs
            empty = 0; tot = 0
            for I in CartesianIndices(ntuple(_ -> 33, ctx[:D]))
                y = HCC._wrap_point(pg_before, Space{ctx[:D], T}(ntuple(d -> ev.center[d] - 2Rw + (I[d] - 0.5) * 4Rw / 33, ctx[:D])))
                n = HCC._collect_near!(pg_before, y, Rw^2, idx, dv); tot += 1; n == 0 && (empty += 1)
            end
            println("HOLECHECK step=$(k+1) center=$(Tuple(round.(ev.center, digits=2))) uncovered quadrature points (no particle within R): $empty of $tot")
        end
    end
    if evs !== nothing && maximum(abs, r.mass_residual) > 1e-11 * maximum(abs, m0)
        for ev in evs
            println("EVENT step=$(k+1) residual=$(r.mass_residual[1]) center=$(Tuple(ev.center)) removed=$(length(ev.removed)) added=$(length(ev.added))")
        end
    end
    # How much of Σ V_i U_i is lost when the estimator replaces the remap's volumes (measured on a copy of the
    # volumes, then restored, so the run itself is unchanged)
    if r.inserted + r.removed > 0 && ts.volume_estimator !== nothing
        Vkeep = copy(collect(view(pg.core.volumes, 1:pg.meta.N))); mr = mass(pg)
        HyperCloud.update_volumes!(pg, ts.volume_estimator)
        global sum_revol += mass(pg) - mr
        if get(ENV, "ADAPT_KEEP_CONTENT", "0") == "1"
            # content-preserving variant: keep m_i = V_i U_i of the remap with the re-estimated volumes
            Vn = view(pg.core.volumes, 1:pg.meta.N); U = view(pg.rhos, 1:pg.meta.N)
            for i in eachindex(Vn); U[i] = U[i] * (Vkeep[i] / Vn[i]); end
        else
            copyto!(view(pg.core.volumes, 1:pg.meta.N), Vkeep)
        end
    end
    global t += dt; global k += 1
    global ins += r.inserted; global rem += r.removed
    global sum_remap += r.mass_residual; global sum_step += ma - mb
    global max_remap = max(max_remap, maximum(abs, r.mass_residual))
    println(io, k, ",", t, ",", pg.meta.N, ",", r.inserted, ",", r.removed, ",", r.mass_residual[1], ",", (ma - mb)[1], ",", vol(pg), ",", mass(pg)[1])
    if t >= next_snap - 1e-12
        push!(snaps, (Float64(t), collect(view(pg.core.positions, 1:pg.meta.N)), snapstate(pg)))
        global next_snap += tmax / nsnap
    end
    any(x -> !isfinite(x), mass(pg)) && (global status = "non-finite"; break)
end
close(io)
U = collect(view(pg.rhos, 1:pg.meta.N))
IS_BGK && (U = [HyperCloud.moments(eq, u) for u in U])
minu1 = minimum(u -> u[1], U)
summary = Dict("name" => name, "method" => method, "tag" => tag, "params" => p, "status" => status, "t_end" => Float64(t),
               "steps" => k, "N0" => N0, "N_end" => pg.meta.N, "inserted" => ins, "removed" => rem,
               "mass0" => collect(m0), "mass_end" => collect(mass(pg)), "sum_remap" => collect(sum_remap), "sum_step" => collect(sum_step),
               "max_abs_remap" => Float64(max_remap), "min_u1" => Float64(minu1), "runtime" => elapsed, "snapshots" => snaps)
JLD2.save(joinpath(out, tag * ".jld2"), summary)
rel(x) = x ./ max.(abs.(collect(m0)), 1e-300)
@printf("RESULT tag=%s status=%s t=%.4g/%.4g steps=%d N=%d->%d inserted=%d removed=%d min_u1=%.4e runtime=%.1fs\n",
        tag, status, t, tmax, k, N0, pg.meta.N, ins, rem, minu1, elapsed)
println("RESULT tag=$tag rel_mass_change=", join((@sprintf("%+.3e", x) for x in rel(collect(mass(pg) - m0))), ","),
        " rel_remap_sum=", join((@sprintf("%+.3e", x) for x in rel(collect(sum_remap))), ","),
        " rel_step_sum=", join((@sprintf("%+.3e", x) for x in rel(collect(sum_step))), ","),
        @sprintf(" max_abs_remap=%.3e", max_remap), " rel_revolume_sum=", join((@sprintf("%+.3e", x) for x in rel(collect(sum_revol))), ","))

# Quadrature check (linear advection only, exact solution u(x, t) = u0(x - a t), periodic): the discrete mass of the
# EXACT solution on the final cloud with the final volumes, Σ V_i u_exact(x_i), against the exact mass. If this alone
# deviates, the "mass change" above is a property of the volume estimate on the moved cloud, not of the scheme.
if p[:PDE_name] == :linear
    ic = HyperCloud.build_ic(ns(:IC), ctx)
    a = p[:PDE_velocities][1][1]; lo, hi = p[:Grid_mins][1], p[:Grid_maxs][1]; L = hi - lo
    X = collect(view(pg.core.positions, 1:pg.meta.N)); V = collect(view(pg.core.volumes, 1:pg.meta.N))
    ue(x) = ic(Space{1, T}(lo + mod(x[1] - a * t - lo, L)))[1]
    qe = sum(V[i] * ue(X[i]) for i in eachindex(X))
    @printf("RESULT tag=%s quad_exact_rel=%+.3e scheme_vs_exact_rel=%+.3e (relative to the exact mass %.6f)\n",
            tag, (qe - m0[1]) / m0[1], (mass(pg)[1] - qe) / m0[1], m0[1])
end

# Smooth-wave error (experiments/BGK/bgk_wave.jl and any IC with a `wave_density(x)` and velocity `WAVE_U`): the
# Euler-limit solution is the initial density translated with u, so L1(ρ) = Σ V_i |ρ_i − ρ0(x_i − u t)|.
if isdefined(Main, :wave_density)
    lo, hi = p[:Grid_mins][1], p[:Grid_maxs][1]; L = hi - lo
    X = collect(view(pg.core.positions, 1:pg.meta.N)); V = collect(view(pg.core.volumes, 1:pg.meta.N))
    werr = sum(V[i] * abs(U[i][1] - Main.wave_density(lo + mod(X[i][1] - Main.WAVE_U * t - lo, L))) for i in eachindex(X))
    @printf("RESULT tag=%s wave_L1_rho=%.4e\n", tag, werr)
end

# Interior-only content (non-periodic problems: the ghosts hold boundary states and are excluded); for Sod-type
# problems the boundary fluxes vanish until the waves reach the ends, so this is the conservation defect.
let n = pg.meta.N, ib = collect(view(pg.core.is_boundary, 1:n)), V = collect(view(pg.core.volumes, 1:n)), F = collect(view(pg.rhos, 1:n))
    mi = sum(V[i] * F[i] for i in 1:n if !ib[i])
    mi = IS_BGK ? HyperCloud.moments(eq, mi) : mi
    println("RESULT tag=$tag interior_content=", join((@sprintf("%.10e", x) for x in mi), ","), " interior_content0=", join((@sprintf("%.10e", x) for x in mint0), ","),
            " rel_interior_change=", join((@sprintf("%+.3e", (mi[c] - mint0[c]) / max(abs(mint0[c]), abs(mint0[1]))) for c in eachindex(mi)), ","))
end

# Vortex error (IsentropicVortex / bgk_vortex, one period: the exact solution is the initial condition)
if isdefined(Main, :vortex_state) && p[:tmax] == 10.0 && p[:PDE_D] == 2
    X = collect(view(pg.core.positions, 1:pg.meta.N)); V = collect(view(pg.core.volumes, 1:pg.meta.N))
    verr = sum(V[i] * abs(U[i][1] - Main.vortex_state(X[i][1], X[i][2])[1]) for i in eachindex(X))
    qex = sum(V[i] * Main.vortex_state(X[i][1], X[i][2])[1] for i in eachindex(X))
    @printf("RESULT tag=%s vortex_L1_rho=%.4e quad_exact_rel=%+.3e scheme_vs_exact_rel=%+.3e\n", tag, verr, (qex - m0[1]) / m0[1], (mass(pg)[1] - qex) / m0[1])
end
