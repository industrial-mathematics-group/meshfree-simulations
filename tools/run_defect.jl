# Runs the conservation-defect experiments (experiments/ConservationDefect) through the PDEStudioCore pipeline
# (simulations, statistics, post-processing that attaches :conservation_defect) and prints one markdown table per
# experiment from the stored statistics. CPU, the configs' own precision (Float64).
#
#     julia -t 12 --project=env_headless tools/run_defect.jl defect_floor defect_shock ... [> tables.md]
#
# Environment variables:
#   DEFECT_SAVE     = directory of the result files (default results/defect)
#   DEFECT_SKIP_RUN = 1 only tabulates the stored runs (missing ones are reported)
#   DEFECT_SERIAL   = 1 runs the tasks one after the other (default: in parallel over the Julia threads,
#                     each run serial; the 1-D runs are too small to profit from threaded particle loops)
#   DEFECT_TABLES   = file the tables are also written to
using Logging, Printf
using HyperCloud, StaticArrays

global_logger(ConsoleLogger(stderr, Logging.Warn))
set_target_module!(Main)
set_save_path!(get(ENV, "DEFECT_SAVE", joinpath(@__DIR__, "..", "results", "defect")))
set_stat_preset!("hyperbolic")
HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")
enable_progress_bar!(false)
set_backend!(PolyesterBackend())

const PARALLEL = get(ENV, "DEFECT_SERIAL", "0") != "1"

function load_config(name)
    path = HyperCloud.resolve_config_path(Main, name)
    Base.include(Main, path)
    return Base.invokelatest(Main.get_config)
end

# Parameters that tell the runs of one config apart, besides the method and N
const SWEPT = (:Grid_randomness_factor, :Weight_range, :Weight_alpha)
setting(p) = join(("$(k == :Grid_randomness_factor ? "θ" : k == :Weight_range ? "κ" : "α") = $(first(p[k]))" for k in SWEPT), ", ")

function method_label(p)
    variant = get(p, :MOOD_criterion, :none) != :none ? p[:MOOD_strategy] : get(p, :Limiter_name, :none)
    return review_label(p[:Scheme_order], variant, p[:Flux_name], p[:Scheme_MLS_order])
end

# Least-squares slope of −log|d| over log N (the observed order of the defect) for N ≥ Nmin
function fitted_rate(Ns, ds; Nmin = 100)
    sel = [i for i in eachindex(Ns) if Ns[i] >= Nmin && isfinite(ds[i]) && ds[i] != 0]
    length(sel) < 3 && return NaN
    x = log.(Ns[sel]); y = -log.(abs.(ds[sel]))
    return sum((x .- sum(x) / length(x)) .* (y .- sum(y) / length(y))) / sum(abs2, x .- sum(x) / length(x))
end

fmt(x) = isnan(x) ? "–" : @sprintf("%+.2e", x)

function table(name, config)
    rows = Dict{Tuple{String, String}, Vector{Tuple{Int, Float64, Float64}}}()   # (method, setting) => (N, defect, floor)
    missing_runs = 0
    for p in generate_tasks(config)
        key = (method_label(p), setting(p))
        N = p[:Grid_Ns][1]
        d = NaN; fl = NaN
        if does_sim_data_exist(p, Val(:raw))
            data = load_sim_data(p, :raw)
            st = data.stats
            if haskey(st, :final_conservation_defect) && data.t[end] >= p[:tmax] * (1 - 1e-6)
                d = st[:final_conservation_defect][1]
                h = prod(data.domain.spacing[1:end-1])
                m = sum(data.u[end])[1] * h
                fl = haskey(p, :IC_uL) ? -h * p[:IC_uL][1] / (2m) : NaN      # half a particle of the inflow state
            end
        else
            missing_runs += 1
        end
        push!(get!(rows, key, Tuple{Int, Float64, Float64}[]), (N, d, fl))
    end
    io = IOBuffer()
    println(io, "\n## $name\n")
    missing_runs > 0 && println(io, "($missing_runs runs missing)\n")
    show_N = (100, 316, 1000, 3163)
    D = config.shared[:PDE_D]
    Nshow = D == 1 ? show_N : (32, 63, 100, 200)
    println(io, "| Method | setting | ", join(("N = $N" for N in Nshow), " | "), " | max \\|d\\| for N ≥ $(D == 1 ? 631 : 100) | fitted rate (N ≥ 100) | −h u_L/(2m) at N_max |")
    println(io, "|", repeat("---|", length(Nshow) + 5))
    for key in sort(collect(keys(rows)))
        r = sort(rows[key]); Ns = first.(r); ds = [x[2] for x in r]
        at(N) = (i = findfirst(==(N), Ns); i === nothing ? NaN : ds[i])
        big = [abs(ds[i]) for i in eachindex(Ns) if Ns[i] >= (D == 1 ? 631 : 100)]
        failed = count(isnan, ds)
        println(io, "| $(key[1]) | $(key[2]) | ", join((fmt(at(N)) for N in Nshow), " | "), " | ",
                isempty(big) ? "–" : @sprintf("%.2e", maximum(big)), failed > 0 ? " ($failed failed)" : "", " | ",
                isnan(fitted_rate(Ns, ds)) ? "–" : @sprintf("%.2f", fitted_rate(Ns, ds)), " | ", fmt(r[end][3]), " |")
    end
    return String(take!(io))
end

out = IOBuffer()
for name in ARGS
    config = load_config(name)
    n = length(generate_tasks(config))
    if get(ENV, "DEFECT_SKIP_RUN", "0") != "1"
        @warn "Running $name: $n runs (parallel: $PARALLEL)"
        t0 = time()
        Base.invokelatest(() -> run_all_simulations(config; calculate_stats = true, parallel = PARALLEL))
        @warn @sprintf("%s finished in %.1f min", name, (time() - t0) / 60)
    end
    t = Base.invokelatest(table, name, config)
    print(t); flush(stdout)
    print(out, t)
end
haskey(ENV, "DEFECT_TABLES") && write(ENV["DEFECT_TABLES"], String(take!(out)))
