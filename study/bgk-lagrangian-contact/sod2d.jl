# Planar 2-D BGK Sod (experiments/BGK/bgk_sod2d.jl): the density averaged over the particles in x-bins of width 1/N_x,
# then the measures of the 1-D runs against the exact Euler solution at t = 0.2: L1 of the binned profile, plateau errors
# (right [0.21, 0.33], left [0.01, 0.16]), contact width (0.29 < ρ < 0.40 in [-0.02, 0.33]) and the shock position.
#     julia --project=env_web study/bgk-lagrangian-contact/sod2d.jl SAVE_DIR NX NY METHOD...
include(joinpath(@__DIR__, "common.jl"))
global_logger(ConsoleLogger(stderr, Logging.Error))

function sod2d_metrics(d, nx)
    xc, ρ = binned_profile(d, nx)
    ok = .!isnan.(ρ); xc = xc[ok]; ρ = ρ[ok]; Δx = 1 / nx
    L1 = sum(abs.(ρ .- ρexact.(xc, d.t[end]))) * Δx
    sel(a, b) = findall(x -> a <= x <= b, xc)
    plat = sum(abs.(ρ[sel(0.21, 0.33)] ./ ρsR .- 1)) / length(sel(0.21, 0.33))
    left = sum(abs.(ρ[sel(0.01, 0.16)] ./ ρsL .- 1)) / length(sel(0.01, 0.16))
    width = count(i -> 0.29 < ρ[i] < 0.40, sel(-0.02, 0.33)) * Δx
    return (; L1, plat, left, width, shock = shock_position((; x = xc, ρ)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    save, nx, ny = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
    c = get_study_config("bgk_sod2d")
    println("| Variant | L1 | plateau err | left err | contact width | shock x |")
    println("|---|---|---|---|---|---|")
    for m in ARGS[4:end]
        set_save_path!(isabspath(save) ? save : joinpath(ROOT, save))
        p = HyperCloud.PDEStudioCore.assemble_params(c.shared, c.methods, Symbol(m)); p[:Grid_Ns] = (nx, ny)
        does_sim_data_exist(p, Val(:raw)) || (println("| $m | missing |"); continue)
        r = sod2d_metrics(load_sim_data(p, :raw), nx)
        @printf("| %s | %.2e | %.1f %% | %.1f %% | %.3f | %.4f |\n", m, r.L1, 100r.plat, 100r.left, r.width, r.shock)
    end
end
