# Tables of the 2-D Euler benchmarks (experiments/Benchmarks2D) from the stored runs: per config and method the time
# reached, NaNs, minimum density and pressure, steps and runtime, and the benchmark's own statistic (Noh: L1/L2 error of
# ρ and L1 of ρu against the exact solution (also Gresho); implosion: x ↔ y symmetry error; Sedov: shock radius, predicted
# radius, angular spread, maximum density; Kelvin–Helmholtz: amplitude of the seeded mode). Runs are looked up through the configs (missing runs are listed as such).
#
#     julia --project=env_headless tools/bench_table.jl [--save DIR] NAME... [> tables.md]
using HyperCloud, StaticArrays, Printf
args = copy(ARGS); save = "."
if !isempty(args) && args[1] == "--save"; popfirst!(args); save = popfirst!(args); end
set_target_module!(Main); set_save_path!(abspath(save)); HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")

function load_config(name)
    c = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
    c isa SimulationConfig || (c = Base.invokelatest(Main.get_config))
    return c
end
label(t) = t[:Scheme_order] == 1 ? "RK2MUSCL1" :
    Base.invokelatest(Main.review_label, t[:Scheme_order], get(t, :Limiter_name, :none) == :VK ? :VK :
                      get(t, :MOOD_criterion, :none) == :none ? :none : t[:MOOD_strategy], t[:Flux_name], t[:Scheme_MLS_order])
fmt(x) = @sprintf("%.3e", x)

for name in args
    c = load_config(name)
    println("\n## $name\n")
    extra = ""
    rows = String[]
    for t in Base.invokelatest(generate_tasks, c)
        lab = label(t)
        if !does_sim_data_exist(t, Val(:raw))
            push!(rows, "| $lab | no data (run failed or not fetched) | | | | | | |"); continue
        end
        d = load_sim_data(t, :raw); st = d.stats; u = d.u[end]; γ = t[:PDE_gamma]
        pr(s) = (γ - 1) * (s[4] - (s[2]^2 + s[3]^2) / (2s[1]))
        nan = count(s -> any(isnan, s), u)
        ex = if haskey(st, :sedov_shock)
            s = st[:sedov_shock]; extra = "R_shock / R_pred / spread / ρ_max"; @sprintf("%.3f / %.3f / %.3f / %.2f", s...)
        elseif haskey(st, :final_kh_mode)
            extra = "KH mode amplitude M(t_end) (M(0) = 0.01)"; fmt(st[:final_kh_mode][1])
        elseif haskey(st, :final_symmetry_error)
            extra = "symmetry error (final)"; fmt(st[:final_symmetry_error][1])
        elseif c.reference !== nothing && haskey(st, :l1error)
            extra = "L1(ρ) / L2(ρ) / L1(ρu) vs exact"
            fmt(last(st[:l1error])[1]) * " / " * fmt(last(st[:l2error])[1]) * " / " * fmt(last(st[:l1error])[2])
        else
            ""
        end
        push!(rows, @sprintf("| %s | %.4g / %.4g | %d | %.4f | %.3e | %d | %.0f | %s |", lab, d.t[end], t[:tmax], nan,
                             minimum(s -> s[1], u), minimum(pr, u), Int(get(st, :k_step, [0])[1]), get(st, :Runtime, [NaN])[1], ex))
    end
    println("| Method | t / t_end | NaN | min ρ | min p | steps | runtime [s] | $(extra) |")
    println("|---|---|---|---|---|---|---|---|")
    foreach(println, rows)
end
