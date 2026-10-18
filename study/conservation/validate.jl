# The standalone :gfd solver reproduces the library right-hand side.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
for D in (1, 2), (order, method) in ((1, (:RK2, :MUSCL, 1)), (2, (:RK2, :MUSCL, 2))), pde in (:burgers, :advection)
    eq, pg, ts, _ = setup(D = D, N = D == 1 ? 100 : 30, rf = 0.2, method = method, pde = pde)
    u = [pg.rhos[i][1] for i in 1:pg.meta.N]
    P = Pairs(pg; deg = 1)
    K = similar(u); rhs!(K, P, nothing, eq, u; order = order)
    Kl = library_rhs(ts, eq, pg, u)
    @printf("D=%d order=%d %-9s  max|K - K_lib| / max|K_lib| = %.2e\n", D, order, pde, maximum(abs, K .- Kl) / maximum(abs, Kl))
end
