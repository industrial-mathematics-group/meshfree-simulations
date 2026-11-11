# Figures of the study (CairoMakie): density profiles, the factor overview, shock position and mass, MOOD activity.
#     julia --project=env_web study/bgk-lagrangian-contact/figures.jl
# Needs the runs in data/contact and data/contact_alepatch and the MOOD profiles (mood_profile.jl), see README.md.
include(joinpath(@__DIR__, "common.jl"))
using PDEStudio.CairoMakie
global_logger(ConsoleLogger(stderr, Logging.Error))
const IMG = joinpath(@__DIR__, "img"); mkpath(IMG)
const C = (blue = colorant"#2a78d6", orange = colorant"#eb6834", aqua = colorant"#1baf7a", yellow = colorant"#eda100",
           exact = colorant"#52514e", grid = colorant"#e4e3df")
set_theme!(Theme(fontsize = 15, Axis = (xgridcolor = C.grid, ygridcolor = C.grid, topspinevisible = false,
                                        rightspinevisible = false)))
const BGK, PATCH = "data/contact", "data/contact_alepatch"
const SERIES = [("EPD2 static", BGK, "EPD2", C.blue, "EPD2 static"),
                ("EPD2", BGK, "EPD2", C.orange, "EPD2 Lagrangian"),
                ("EPD2", PATCH, "EPD2", C.aqua, "EPD2 Lagrangian + ALE flux patch"),
                ("VK", BGK, "VK", C.yellow, "VK Lagrangian")]

# 1. density profiles between rarefaction and shock
let fig = Figure(size = (1200, 480))
    for (j, N) in enumerate((800, 3200))
        ax = Axis(fig[1, j], title = "BGK Sod, τ = 1e-5, N = $N, t = 0.2", xlabel = "x", ylabel = j == 1 ? "ρ" : "",
                  limits = ((0.0, 0.42), (0.1, 0.5)))
        xe = range(0.0, 0.42, length = 2000)
        lines!(ax, xe, ρexact.(xe, 0.2), color = C.exact, linestyle = :dash, linewidth = 1.5, label = "exact Euler")
        for (m, save, _, col, lab) in SERIES
            s = snapshot(load_run("bgk_sod_contact", save, m, N))
            lines!(ax, s.x, s.ρ, color = col, linewidth = 2, label = lab)
        end
        j == 2 && axislegend(ax, position = :rt, framevisible = false)
    end
    save(joinpath(IMG, "profiles.png"), fig)
end

# 2. factor overview: plateau error at N = 800 per variant (metrics CSVs of metrics.jl)
function read_metrics(path)
    rows = Dict{Tuple{String, Int}, NamedTuple}()
    for l in Iterators.drop(eachline(path), 1)
        f = split(l, ','); isempty(f[3]) && continue
        rows[(f[1], parse(Int, f[2]))] = (L1 = parse(Float64, f[3]), plat = parse(Float64, f[4]), width = parse(Float64, f[6]))
    end
    return rows
end
let base = read_metrics(joinpath(ROOT, BGK, "metrics_bgk_sod_contact.csv")),
    patch = read_metrics(joinpath(ROOT, PATCH, "metrics_bgk_sod_contact.csv"))
    names = [k[1] for k in keys(base) if k[2] == 800]
    sort!(names, by = n -> base[(n, 800)].plat)
    fig = Figure(size = (900, 760))
    ax = Axis(fig[1, 1], title = "Contact: mean plateau error between contact and shock, N = 800",
              xlabel = "mean |ρ/ρ* − 1| for x ∈ [0.21, 0.33]  (%)", yticks = (1:length(names), names), xscale = log10,
              xticks = [0.1, 0.3, 1, 3, 10], xtickformat = v -> string.(v))
    vlines!(ax, [100base[("EPD2 static", 800)].plat], color = C.blue, linestyle = :dot, linewidth = 1.5)
    scatter!(ax, [100base[(n, 800)].plat for n in names], 1:length(names), color = C.orange, markersize = 11,
             label = "unpatched")
    ip = [i for (i, n) in enumerate(names) if haskey(patch, (n, 800))]
    scatter!(ax, [100patch[(names[i], 800)].plat for i in ip], ip, color = C.aqua, marker = :diamond, markersize = 12,
             label = "ALE flux patch")
    axislegend(ax, position = :rb, framevisible = false)
    text!(ax, 100base[("EPD2 static", 800)].plat, length(names) + 0.3, text = "EPD2 static", color = C.exact,
          fontsize = 12, align = (:left, :bottom))
    save(joinpath(IMG, "factors.png"), fig)
end

# 3. shock position over N and mass over time
let fig = Figure(size = (1200, 460))
    Ns = [200, 400, 800, 1600, 3200]
    ax1 = Axis(fig[1, 1], title = "Shock position at t = 0.2", xlabel = "N", ylabel = "x of the shock", xscale = log2,
               xticks = (Ns, string.(Ns)))
    hlines!(ax1, [0.2S], color = C.exact, linestyle = :dash, linewidth = 1.5, label = "exact")
    sh = [("bgk_sod_contact", BGK, "EPD2 static", C.blue, :circle, "BGK EPD2 static"),
          ("bgk_sod_contact", BGK, "EPD2", C.orange, :circle, "BGK EPD2 Lagrangian"),
          ("bgk_sod_contact", PATCH, "EPD2", C.aqua, :circle, "BGK EPD2 Lagrangian + patch"),
          ("euler_sod_contact", BGK, "EPD2", C.yellow, :rect, "Euler EPD2 Lagrangian")]
    for (cfg, save, m, col, mk, lab) in sh
        xs = [shock_position(snapshot(load_run(cfg, save, m, N))) for N in Ns]
        scatterlines!(ax1, Ns, xs, color = col, marker = mk, markersize = 10, linewidth = 2, label = lab)
    end
    axislegend(ax1, position = :rt, framevisible = false, labelsize = 12)
    ax2 = Axis(fig[1, 2], title = "Mass Σ hᵢ ρᵢ, N = 3200", xlabel = "t", ylabel = "relative change from t = 0")
    for (cfg, save, m, col, mk, lab) in sh
        d = load_run(cfg, save, m, 3200)
        m0 = sum(snapshot(d, 1).h .* snapshot(d, 1).ρ)
        ms = [sum(snapshot(d, k).h .* snapshot(d, k).ρ) / m0 - 1 for k in eachindex(d.t)]
        lines!(ax2, d.t, ms, color = col, linewidth = 2, label = lab)
    end
    axislegend(ax2, position = :lt, framevisible = false, labelsize = 12)
    save(joinpath(IMG, "shock_mass.png"), fig)
end

# 4. MOOD activity along x (fraction of steps with the EPD2 edge order below K, Eulerian bins)
let fig = Figure(size = (900, 420))
    ax = Axis(fig[1, 1], title = "MOOD activity, BGK EPD2, N = 800, t ∈ [0, 0.2]", xlabel = "x",
              ylabel = "fraction of steps with edge order < K", limits = ((-0.35, 0.45), (0, nothing)))
    for (tag, col, lab) in (("epd2s", C.blue, "static"), ("epd2", C.orange, "Lagrangian"),
                            ("epd2p", C.aqua, "Lagrangian + ALE flux patch"))
        rows = [parse.(Float64, split(l, ',')) for l in Iterators.drop(eachline(joinpath(ROOT, BGK, "mood_$(tag)_bins.csv")), 1)]
        lines!(ax, first.(rows), getindex.(rows, 3), color = col, linewidth = 2, label = lab)
    end
    vlines!(ax, [0.2ust, 0.2S], color = C.exact, linestyle = :dash, linewidth = 1)
    text!(ax, [0.2ust, 0.2S], [0, 0], text = ["contact", "shock"], color = C.exact, fontsize = 12,
          align = (:right, :bottom), offset = (-3, 3))
    axislegend(ax, position = :lt, framevisible = false)
    save(joinpath(IMG, "mood_activity.png"), fig)
end
println("figures written to $IMG")
