# Contact metrics of the 1-D Sod runs of a config (bgk_sod_contact, euler_sod_contact) at t = 0.2, per method and N:
#   L1     the :l1error statistic (density, against the exact Euler solution) at the final time
#   plat   mean |ρ/ρ* − 1| between contact and shock (x ∈ [0.21, 0.33], ρ* = 0.26557)
#   left   mean |ρ/ρ*_L − 1| between rarefaction tail and contact (x ∈ [0.01, 0.16], ρ*_L = 0.42632)
#   width  length of x ∈ [−0.02, 0.33] where 0.29 < ρ < 0.40 (the smeared contact; sum of particle spacings)
# Exact Sod at t = 0.2: rarefaction tail x = −0.0141, contact x = 0.1855, shock x = 0.3504.
#     julia --project=env_web study/bgk-lagrangian-contact/metrics.jl CONFIG SAVE_DIR [N ...]
# Prints a markdown table (N = the largest given N, or 800 and 3200 with the L1 order) and writes SAVE_DIR/metrics_CONFIG.csv.
using Logging, Printf
using HyperCloud, StaticArrays
global_logger(ConsoleLogger(stderr, Logging.Warn))
const ROOT = abspath(joinpath(@__DIR__, "..", ".."))
name, save = ARGS[1], ARGS[2]
Ns = length(ARGS) > 2 ? parse.(Int, ARGS[3:end]) : [200, 400, 800, 1600, 3200]
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = ROOT
set_save_path!(isabspath(save) ? save : joinpath(ROOT, save))
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))

function metrics(d)
    x = [p[1] for p in d.x[end]]; ρ = [u[1] for u in d.u[end]]
    perm = sortperm(x); x = x[perm]; ρ = ρ[perm]
    h = [i == 1 ? x[2] - x[1] : i == length(x) ? x[end] - x[end - 1] : (x[i + 1] - x[i - 1]) / 2 for i in eachindex(x)]
    sel(a, b) = findall(xi -> a <= xi <= b, x)
    ip, il = sel(0.21, 0.33), sel(0.01, 0.16)
    plat = sum(abs.(ρ[ip] ./ 0.26557 .- 1)) / length(ip)
    left = sum(abs.(ρ[il] ./ 0.42632 .- 1)) / length(il)
    iw = [i for i in sel(-0.02, 0.33) if 0.29 < ρ[i] < 0.40]
    width = sum(h[iw]; init = 0.0)
    r27 = sum(ρ[sel(0.26, 0.28)]) / length(sel(0.26, 0.28))
    return (; plat, left, width, r27)
end

rows = []
for m in config.defaults, N in Ns
    p = HyperCloud.PDEStudioCore.assemble_params(config.shared, config.methods, m)
    p[:Grid_Ns] = (N,)
    if !does_sim_data_exist(p, Val(:raw))
        push!(rows, (; method = string(m), N, ok = false)); continue
    end
    d = load_sim_data(p, :raw)
    L1 = haskey(d.stats, :l1error) ? last(d.stats[:l1error])[1] : NaN
    push!(rows, (; method = string(m), N, ok = true, L1, metrics(d)..., steps = d.stats[:k_step][1]))
end
open(joinpath(isabspath(save) ? save : joinpath(ROOT, save), "metrics_$name.csv"), "w") do io
    println(io, "method,N,L1,plat,left,width,rho027,steps")
    for r in rows
        r.ok || (println(io, "$(r.method),$(r.N),,,,,,"); continue)
        @printf(io, "%s,%d,%.6e,%.6e,%.6e,%.6e,%.6f,%d\n", r.method, r.N, r.L1, r.plat, r.left, r.width, r.r27, r.steps)
    end
end
missing_runs = [r for r in rows if !r.ok]
isempty(missing_runs) || println("missing: ", join(("$(r.method) N=$(r.N)" for r in missing_runs), ", "))
Nhi = maximum(Ns); Nlo = Nhi ÷ 4
println("| Variant | L1 N=$Nlo | L1 N=$Nhi | order | plateau err N=$Nlo | plateau err N=$Nhi | left err N=$Nhi | ρ(0.27) N=$Nhi | contact width N=$Nlo / $Nhi |")
println("|---|---|---|---|---|---|---|---|---|")
for m in config.defaults
    lo = findfirst(r -> r.method == string(m) && r.N == Nlo && r.ok, rows)
    hi = findfirst(r -> r.method == string(m) && r.N == Nhi && r.ok, rows)
    (lo === nothing || hi === nothing) && continue
    a, b = rows[lo], rows[hi]
    @printf("| %s | %.2e | %.2e | %.2f | %.1f %% | %.1f %% | %.1f %% | %.4f | %.4f / %.4f |\n", m, a.L1, b.L1,
            log(a.L1 / b.L1) / log(Nhi / Nlo), 100a.plat, 100b.plat, 100b.left, b.r27, a.width, b.width)
end
