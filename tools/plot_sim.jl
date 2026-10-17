# Final-time density (component 1) of every run of an experiment as a scatter image, using the same
# SIM_* selection as run_sim.jl:  julia --project=env_plot tools/plot_sim.jl NAME [snapshot index]
# Writes figures/<NAME>_<label>.png.
using HyperCloud, StaticArrays, CairoMakie, Printf
set_target_module!(Main)
set_save_path!(get(ENV, "SIM_SAVE", joinpath(@__DIR__, "..", "results")))
HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")
include(joinpath(@__DIR__, "config.jl"))
mkpath(joinpath(@__DIR__, "..", "figures"))
comp = get(ENV, "PLOT_COMPONENT", "1")     # a component index or "p" (Euler pressure)
value(s, task) = comp == "p" ? (task[:PDE_gamma] - 1) * (s[end] - sum(abs2, s[2:end-1]) / (2s[1])) : s[parse(Int, comp)]
for task in generate_tasks(config)
    does_sim_data_exist(task, Val(:raw)) || continue
    d = load_sim_data(task, :raw)
    n = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : length(d.t)
    x, u = d.x[n], d.u[n]
    v = [value(s, task) for s in u]
    imin = argmin(v); println("min $(comp) = $(v[imin]) at $(x[imin]); ", count(<(0), v), " negative")
    lab = label(task)
    D = length(first(x))
    fig = if D == 1
        f = Figure(size = (900, 500)); ax = Axis(f[1, 1], title = "$NAME  $lab  t=$(round(d.t[n], digits = 4))")
        o = sortperm(first.(x)); lines!(ax, first.(x)[o], v[o]); f
    else
        xs, ys = first.(x), last.(x)
        aspect = (maximum(xs) - minimum(xs)) / (maximum(ys) - minimum(ys))
        f = Figure(size = (round(Int, 900 * min(aspect, 3)) + 150, 900))
        ax = Axis(f[1, 1], title = "$NAME  $lab  t=$(round(d.t[n], digits = 4))", aspect = DataAspect())
        ms = 1.2 * 900 / sqrt(length(x) / aspect)
        sc = scatter!(ax, xs, ys, color = v, markersize = ms, colormap = :turbo,
                      colorrange = extrema(filter(isfinite, v)))
        Colorbar(f[1, 2], sc); f
    end
    file = joinpath(@__DIR__, "..", "figures", replace("$(NAME)_$(lab)_T$(task[:T])_$(comp)$(get(ENV, "PLOT_TAG", "")).png", r"[^A-Za-z0-9_.()=,x-]" => "_"))
    save(file, fig)
    println("wrote $file")
end
