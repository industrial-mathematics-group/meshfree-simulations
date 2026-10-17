include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
for (pde, method, R) in ((:burgers, (:RK2, :MUSCL, 2), 2.5), (:burgers, (:RK2, :MUSCL, 1), 2.5), (:advection, (:RK2, :MUSCL, 1), 2.5), (:burgers, (:RK4, :MUSCL, 4), 5.5)), rf in (0.0, 0.2)
    t = @elapsed begin
        eq, pg, ts, _ = setup(D = 2, N = 64, rf = rf, method = method, pde = pde, R = R)
        r = evolve!(eq, pg, ts; tmax = 0.5)
    end
    @printf("%-9s %-10s R=%.1f rf=%.1f  dVU=%+.2e dU=%+.2e steps=%d  [%.3f, %.3f]  %.1fs\n", pde, join(string.(method)), R, rf, r.dVU, r.dU, r.steps, r.umin, r.umax, t)
end
