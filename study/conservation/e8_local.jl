# E8: fully local conservative variants built from intercept coefficients (no global solve):
#   gfd_i : -2 Σ c_ij·(F*_ij - F_i) with intercept c (Σ_j c_ij = 0)
#   sym_i : -(1/V_i) Σ (V_i c_ij - V_j c_ji)·F*_ij   (exactly conservative, normal Riemann solver)
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf, LinearAlgebra
u0(D) = D == 1 ? (x -> 1 + 0.5 * sinpi(2x[1])) : (x -> 1 + 0.5 * sinpi(2x[1]) * sinpi(2x[2]))
INTERCEPT[] = true; SELF[] = false
for D in (1, 2), R in (2.5, 5.5), rf in (0.05, 0.2), (form, rs) in ((:gfd, :axis), (:sym, :normal))
    for N in (D == 1 ? (100, 400, 1600) : (32, 64, 128))
        eq, pg, ts, _ = setup(D = D, N = N, rf = rf, R = R)
        n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]; V = collect(pg.core.volumes[1:n])
        P = Pairs(pg; deg = 1); A = face_coefficients(P, form)
        col = A === nothing ? 0.0 : maximum(i -> norm(sum(A[k] for k in P.ranges[i])) / sum(norm(A[k]) for k in P.ranges[i]), 1:n)
        ui = [u0(D)(xi) for xi in x]
        us, _ = solve(P, A, eq, pg, ui, 0.2; order = 2, riemann = rs)
        es = sum(V .* abs.(us .- [burgers_exact(u0(D), xi, 0.2) for xi in x]))
        up, _ = solve(P, A, eq, pg, ui, D == 1 ? 1.0 : 0.5; order = 2, limiter = :bj, riemann = rs)
        @printf("D=%d R=%.1f rf=%.2f %-4s N=%4d  max|ΣA|/Σ|A|=%.1e | smooth L1=%.2e | post-shock |dVU|=%.2e\n",
                D, R, rf, form, N, col, es, abs(sum(V .* up) - sum(V .* ui)) / sum(V .* ui))
    end
end
