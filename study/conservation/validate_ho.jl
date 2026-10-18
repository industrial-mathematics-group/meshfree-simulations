include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
for D in (1, 2), o in (2, 3, 4)
    eq, pg, ts, _ = setup(D = D, N = D == 1 ? 100 : 30, rf = 0.2, R = 5.5, method = (:RK4, :MUSCL, o))
    u = [pg.rhos[i][1] for i in 1:pg.meta.N]
    P = PairsHO(pg; q = o - 1); C = gfd_coefficients(P)
    K = similar(u); rhs_ho!(K, P, nothing, C, eq, u)
    Kl = library_rhs(ts, eq, pg, u)
    @printf("D=%d MUSCL%d  max|K - K_lib| / max|K_lib| = %.2e\n", D, o, maximum(abs, K .- Kl) / maximum(abs, Kl))
end
