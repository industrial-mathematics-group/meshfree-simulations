# E5: reconstruction degree vs divergence-fit degree (smooth advection, regular and jittered grid).
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
u0 = x -> 1 + 0.5 * sinpi(2x[1])
for rf in (0.0, 0.2), (rq, dq) in ((1, 1), (1, 2), (2, 2), (2, 3), (3, 2), (3, 3), (4, 4), (4, 5), (2, 1), (1, 3))
    errs = Float64[]; Ns = (50, 100, 200, 400)
    for N in Ns
        eq, pg, ts, _ = setup(D = 1, N = N, rf = rf, R = 5.5, pde = :advection)
        n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]; V = collect(pg.core.volumes[1:n])
        P = PairsHO(pg; q = rq); C = gfd_coefficients(PairsHO(pg; q = dq))
        u, _ = solve_ho(P, nothing, eq, pg, [u0(xi) for xi in x], 0.5; cfl = 0.1, C = C)
        push!(errs, sum(V .* abs.(u .- [u0(xi .- 0.5) for xi in x])))
    end
    eoc = [log(errs[i-1] / errs[i]) / log(Ns[i] / Ns[i-1]) for i in 2:length(Ns)]
    @printf("rf=%.1f reconstruction degree %d, divergence degree %d:  L1 %s   EOC %s\n", rf, rq, dq, join((@sprintf("%.2e", e) for e in errs), " "), join((@sprintf("%.2f", e) for e in eoc), " "))
end
