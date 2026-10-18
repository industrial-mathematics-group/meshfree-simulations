# E7: MLS divergence WITH a constant term (zero row sums) vs the library fit (no constant term).
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf, LinearAlgebra
u0(D) = D == 1 ? (x -> 1 + 0.5 * sinpi(2x[1])) : (x -> 1 + 0.5 * sinpi(2x[1]) * sinpi(2x[2]))
for D in (1, 2), R in (2.5, 5.5), icp in (true,)
    SELF[] = false
    INTERCEPT[] = icp
    for N in (D == 1 ? (100, 400, 1600) : (32, 64, 128))
        eq, pg, ts, _ = setup(D = D, N = N, rf = 0.2, R = R)
        n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]; V = collect(pg.core.volumes[1:n])
        I, J, C = coefficients(pg, 1); a = residual_split(pg, I, J, C, V)
        P = Pairs(pg; deg = 1)                            # Pairs uses `coefficients` → follows INTERCEPT
        ui = [u0(D)(xi) for xi in x]
        us, _ = solve(P, nothing, eq, pg, ui, 0.2; order = 2)
        es = sum(V .* abs.(us .- [burgers_exact(u0(D), xi, 0.2) for xi in x]))
        up, _ = solve(P, nothing, eq, pg, ui, D == 1 ? 1.0 : 0.5; order = 2, limiter = :bj)
        @printf("D=%d R=%.1f intercept=%-5s N=%4d  asym=%.1e rowsum=%.1e | smooth L1=%.2e | post-shock |dVU|=%.2e\n",
                D, R, icp, N, a.asym, a.rowsum, es, abs(sum(V .* up) - sum(V .* ui)) / sum(V .* ui))
    end
end
