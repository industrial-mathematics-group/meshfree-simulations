# E10: library before (env_base) / after (env_study) the divergence changes: post-shock mass drift and
# smooth accuracy. Run with both environments.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
println("== drift after the shock (Burgers, rf = 0.2)")
for D in (1, 2), R in (2.5, 5.5), (method, lim, mood, mls) in (((:RK2, :MUSCL, 1), :none, :none, 0), ((:RK2, :MUSCL, 2), :VK, :none, 0),
                                                             ((:RK4, :MUSCL, 4), :none, :U2, 4), ((:RK4, :MUSCL, 3), :none, :U2, 3))
    N = D == 1 ? 400 : 64
    eq, pg, ts, _ = setup(D = D, N = N, rf = 0.2, R = R, method = method, limiter = lim, mood = mood, mls = mls)
    r = evolve!(eq, pg, ts; tmax = D == 1 ? 1.0 : 0.5)
    @printf("  D=%d R=%.1f %-10s %-4s mood=%-4s MLS=%d  |dVU|=%.2e  [%.3f, %.3f]\n", D, R, join(string.(method)), lim, mood, mls, abs(r.dVU), r.umin, r.umax)
end
println("== smooth advection, L1 at t = 0.5 (rf = 0.2)")
u0(D) = D == 1 ? (x -> 1 + 0.5 * sinpi(2x[1])) : (x -> 1 + 0.5 * sinpi(2x[1]) * sinpi(2x[2]))
for D in (1, 2), R in (2.5, 3.5), o in (2, 3, 4)
    errs = Float64[]; Ns = D == 1 ? (100, 200, 400) : (32, 64, 128)
    for N in Ns
        eq, pg, ts, _ = setup(D = D, N = N, rf = 0.2, R = R, pde = :advection, method = (:RK4, :MUSCL, o))
        evolve!(eq, pg, ts; tmax = 0.5)
        v = D == 1 ? (1.0,) : (1.0, 0.5)
        push!(errs, sum(pg.core.volumes[i] * abs(pg.rhos[i][1] - u0(D)(pg.core.positions[i] .- 0.5 .* v)) for i in 1:pg.meta.N))
    end
    @printf("  D=%d R=%.1f MUSCL%d: %s  EOC %.2f %.2f\n", D, R, o, join((@sprintf("%.2e", e) for e in errs), " "), log2(errs[1] / errs[2]), log2(errs[2] / errs[3]))
end
