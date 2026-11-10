# Figures of the BGK Sod runs (HyperCloudCore 0.8.0 rerun, reports/2026-11-09-bgk-rerun-080.md) with PDEStudio's agent
# API, headless (CairoMakie), into figures/bgk080/. The studio gets a config reconstructed from the stored runs
# (config_from_data.jl), so HyperCloud isn't needed: PDEStudio 0.6 needs PDEStudioCore 0.5, HyperCloud 0.6 still 0.4.
#
# Environment (outside the repo; the repo environments keep PDEStudio 0.5 until HyperCloud supports PDEStudioCore 0.5):
#     julia --project=<scratch> -e 'using Pkg; Pkg.add(["PDEStudio", "HyperCloudCore", "CairoMakie", "StaticArrays"])'
# Data: copies of the 0.8.0 BGK Sod runs, upgraded to the PDEStudioCore 0.5 format, in data/pdestudio05/ (the
# originals in data/ stay readable by the PDEStudioCore 0.4 environments):
#     set_save_path!("data/pdestudio05"); upgrade_sim_data()
# Run:  cd <repo>; julia --project=<scratch> -t 8 tools/studio/bgk_sod_figures.jl
using PDEStudio, HyperCloudCore, StaticArrays, CairoMakie
include(joinpath(@__DIR__, "config_from_data.jl"))
const ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "figures", "bgk080"); mkpath(OUT)
set_save_path!(joinpath(ROOT, "data", "pdestudio05"))
set_stat_preset!("hyperbolic")

lab(m) = string(m[:Time_stepper], "MUSCL", m[:Scheme_order], "(", get(m, :Limiter_name, :none) == :VK ? "VK" : string(m[:MOOD_strategy]),
                ") ", m[:Grid_mover] == :none ? "static" : m[:Adapt_enabled] ? "Lagrangian+adapt" : "Lagrangian")
config = config_from_data(sweep = Dict(:Grid_Ns => :Grid_Ns__1, :PDE_tau => :PDE_tau), label = lab)
@assert check_config(config).missing == 0 "the reconstructed config doesn't find every stored run"

studio_headless!()
launch_plotter(config)

# 1. All nine variants at t = 0.2, N = 800, τ = 1e-5 (particles; lines would join particles in storage order, which
#    particle management doesn't keep sorted), and a zoom on the contact and the shock
studio_set!(:Grid_Ns__1 => 800, :PDE_tau => 1e-5, :t => 0.2, :compare_target => :methods, :compare_columns => 3,
            :plot_style => :particles_1d, :plot_width => 470, :plot_height => 260)
studio_option!(:labels, :x_label, "x"); studio_option!(:labels, :y_label, "ρ")
studio_option!(:plot_style, :marker_size, 3.0)
studio_snapshot(joinpath(OUT, "sod_methods_N800_tau1e-5.png"))
studio_view!(xlims = (0.05, 0.45), ylims = (0.1, 0.5))
studio_snapshot(joinpath(OUT, "sod_methods_N800_tau1e-5_zoom.png"))

# 2. The relaxation time: τ = 1e-5, 1e-3, 1e-2 side by side, three methods overlaid
studio_set!(:methods => Symbol.(["ARS222MUSCL2(VK) static", "ARS222MUSCL2(EPD2) Lagrangian", "ARS233MUSCL3(EPD2) Lagrangian"]))
studio_set!(:compare_target => :PDE_tau, :compare_columns => 3, :plot_style => :particles_lines, :plot_width => 470, :plot_height => 300)
studio_option!(:plot_style, :line_width, 1.2); studio_option!(:plot_style, :marker_size, 4.0)
studio_press!(:update_plot)
studio_snapshot(joinpath(OUT, "sod_tau_N800.png"))

# 3. Final L1(ρ) error against N (CairoMakie: the studio offers statistics only on the x axis, next to Solution)
err = Dict{Tuple{String, Float64}, Vector{Tuple{Int, Float64}}}()
for t in generate_tasks(config)
    push!(get!(err, (lab(t), t[:PDE_tau]), Tuple{Int, Float64}[]), (t[:Grid_Ns][1], last(load_sim_data(t, :raw).stats[:l1error])[1]))
end
fig = Figure(size = (1500, 470), fontsize = 15)
styles = Dict("static" => :solid, "Lagrangian" => :dash, "Lagrangian+adapt" => :dot)
for (j, τ) in enumerate(sort(config.varied[:PDE_tau]))
    ax = Axis(fig[1, j]; xscale = log10, yscale = log10, xlabel = "N", ylabel = j == 1 ? "L1 error of ρ at t = 0.2" : "", title = "τ = $τ")
    for m in string.(config.defaults)
        pts = sort(err[(m, τ)])
        c = occursin("ARS233", m) ? :green : occursin("(VK)", m) ? :black : :blue
        lines!(ax, first.(pts), last.(pts); color = c, linestyle = styles[split(m)[end]], linewidth = 2, label = m)
    end
    j == 1 && lines!(ax, [200, 1600], 0.02 .* ([200, 1600] ./ 200) .^ -1; color = :gray, linestyle = :dashdot, label = "slope −1")
end
Legend(fig[1, 4], fig.content[1]; framevisible = false)
save(joinpath(OUT, "sod_convergence_l1.png"), fig)
println("figures in ", OUT)
