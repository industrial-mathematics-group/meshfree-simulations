# E4: order reduction of odd MUSCL orders: Rusanov dissipation on vs off (central flux), smooth advection.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
u0 = x -> 1 + 0.5 * sinpi(2x[1])
for rf in (0.0, 0.2), o in (2, 3, 4, 5), diss in (1.0, 0.0)
    DISSIPATION[] = diss
    errs = Float64[]; Ns = (50, 100, 200, 400)
    for N in Ns
        eq, pg, ts, _ = setup(D = 1, N = N, rf = rf, R = 5.5, pde = :advection)
        n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]; V = collect(pg.core.volumes[1:n])
        P = PairsHO(pg; q = o - 1)
        u, _ = solve_ho(P, nothing, eq, pg, [u0(xi) for xi in x], 0.5; cfl = 0.1)
        push!(errs, sum(V .* abs.(u .- [u0(xi .- 0.5) for xi in x])))
    end
    eoc = [log(errs[i-1] / errs[i]) / log(Ns[i] / Ns[i-1]) for i in 2:length(Ns)]
    @printf("rf=%.1f MUSCL%d dissipation=%.0f  L1: %s   EOC: %s\n", rf, o, diss, join((@sprintf("%.2e", e) for e in errs), " "), join((@sprintf("%.2f", e) for e in eoc), " "))
end
