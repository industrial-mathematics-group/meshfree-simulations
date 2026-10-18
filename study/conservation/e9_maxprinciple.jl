# E9: first-order maximum principle (box advection / Burgers, 0 ≤ u ≤ 1) of the face form with the
# conservative high-order-constrained coefficients vs the GFD form with a fixed cubic divergence
# (the MOOD fallback situation with MLS_order = 4).
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf, LinearAlgebra, StaticArrays
box(D) = D == 1 ? (x -> 0.3 < x[1] < 0.6 ? 1.0 : 0.0) : (x -> (0.3 < x[1] < 0.6 && 0.3 < x[2] < 0.6) ? 1.0 : 0.0)
# forward Euler steps of the first-order scheme with dt = cfl / max_i (Σ_j |coef_ij| λ / V_i)
function fe_first_order(P, A, C, eq, u, steps, cfl; face)
    n = P.n; K = similar(u)
    for _ in 1:steps
        λ = maximum(abs, u) + 1.0                       # advection speed ≤ 1, Burgers |u| ≤ 1
        stiff = face ? maximum(i -> sum(norm(A[k]) for k in P.ranges[i]) / P.V[i], 1:n) :
                       maximum(i -> 2 * sum(norm(C[k]) for k in P.ranges[i]), 1:n)
        dt = cfl / (λ * stiff)
        rhs_ho!(K, P, face ? A : nothing, C, eq, u; riemann = :normal)   # P.q = 0 → piecewise constant
        u = u .+ dt .* K
    end
    return u
end
for D in (1, 2), pde in (:advection, :burgers), rf in (0.0, 0.2)
    eq, pg, ts, _ = setup(D = D, N = D == 1 ? 200 : 48, rf = rf, R = 5.5, pde = pde)
    n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]
    u0 = [box(D)(xi) for xi in x]
    P3 = PairsHO(pg; q = 3)
    C3 = gfd_coefficients(P3)                           # fixed cubic divergence (MLS_order = 4)
    A3 = face_coefficients_p(P3, 3)                     # conservative, exact to degree 3
    P0 = PairsHO(pg; q = 1); P0 = Harness.PairsHO(P0.n, P0.ranges, P0.J, P0.dist, [zeros(0) for _ in P0.S], P0.rev, P0.V, 0, P0.invL)
    for (name, face, coef) in (("GFD, cubic divergence", false, C3), ("face form, p = 3", true, A3))
        u = fe_first_order(P0, face ? coef : nothing, face ? C3 : coef, eq, copy(u0), 100, 0.9; face)
        @printf("D=%d %-9s rf=%.1f  %-22s  min=%+.2e  max=%.4f  mass change=%+.1e\n", D, pde, rf, name, minimum(u), maximum(u),
                (sum(P0.V .* u) - sum(P0.V .* u0)) / sum(P0.V .* u0))
    end
end
