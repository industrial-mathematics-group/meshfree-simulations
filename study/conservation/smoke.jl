include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
for w in WEIGHTS, rf in (0.0, 0.2)
    eq, pg, ts, _ = setup(D = 1, N = 100, rf = rf, weight = w, alpha = 1.0, R = 2.5)
    r = evolve!(eq, pg, ts; tmax = 1.0)
    @printf("%-12s rf=%.1f  dVU=%+.2e dU=%+.2e steps=%d  [%.3f, %.3f]\n", w, rf, r.dVU, r.dU, r.steps, r.umin, r.umax)
end
