# Position of the BL-quadrant edges at t_end for every 2-D Clain run in the given result directories:
# x where the density along the bottom strip (0.02<y<0.06) first exceeds 0.33, y likewise along the left strip.
using HyperCloud, Printf
for dir in ARGS
    set_save_path!(dir)
    for f in filter(endswith(".jld2"), readdir(joinpath(dir, "data")))
        d = load_sim_data(f; key = "raw")
        p = d.params
        get(p, :PDE_name, nothing) == :euler && get(p, :PDE_D, 0) == 2 && haskey(p, :IC_states) || continue
        x, u = d.x[end], d.u[end]
        bx = sort([(q[1], s[1]) for (q, s) in zip(x, u) if 0.10 < q[2] < 0.15]); ly = sort([(q[2], s[1]) for (q, s) in zip(x, u) if 0.10 < q[1] < 0.15])
        ex = bx[findfirst(t -> t[2] > 0.33, bx)][1]; ey = ly[findfirst(t -> t[2] > 0.33, ly)][1]
        tx = sort([(q[1], s[1]) for (q, s) in zip(x, u) if 0.85 < q[2] < 0.9]); ty = sort([(q[2], s[1]) for (q, s) in zip(x, u) if 0.85 < q[1] < 0.9])
        sx = tx[findfirst(t -> t[2] > 1.0, tx)][1]; sy = ty[findfirst(t -> t[2] > 1.0, ty)][1]
        @printf("%-55s N=%s R=%.1f rf=%s T=%s  x_edge=%.4f y_edge=%.4f  TR shock x=%.4f y=%.4f\n", "$(p[:Time_stepper])$(p[:Scheme_name])$(p[:Scheme_order]) $(get(p,:Limiter_name,:none)) $(get(p,:MOOD_strategy,:-)) CFL=$(p[:Time_CFL])",
                p[:Grid_Ns], p[:Weight_range], p[:Grid_randomness_factor], get(p, :T, Float64), ex, ey, sx, sy)
    end
end
