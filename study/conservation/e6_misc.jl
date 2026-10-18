# E6: (a) weight decay α, (b) MOOD with fixed vs coupled MLS order, (c) time history of the drift,
# (d) the best possible particle weights V for the scheme (least-squares antisymmetry) — library scheme.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf, LinearAlgebra, SparseArrays, StaticArrays, HyperCloud
const HCC = Harness.HCC
part = get(ENV, "PART", "abcd")
if occursin('a', part)
    println("(a) weight decay: 1-D/2-D Burgers, MUSCL2 VK, rf = 0.2, |ΔΣVU| after the shock")
    for D in (1, 2), w in (:exponential, :shifted_exp, :gauss), α in (0.5, 1.0, 2.0, 4.0, 8.0), R in (2.5, 5.5)
        N = D == 1 ? 400 : 64
        eq, pg, ts, _ = setup(D = D, N = N, rf = 0.2, R = R, weight = w, alpha = α, method = (:RK2, :MUSCL, 2), limiter = :VK)
        r = evolve!(eq, pg, ts; tmax = D == 1 ? 1.0 : 0.5)
        @printf("  D=%d %-12s α=%3.1f R=%.1f  |dVU|=%.2e  steps=%d  [%.3f, %.3f]\n", D, w, α, R, abs(r.dVU), r.steps, r.umin, r.umax)
    end
end
if occursin('b', part)
    println("(b) MOOD EPD1 (U2), MUSCL4, R = 5.5: fixed vs coupled MLS order")
    for D in (1, 2), rf in (0.0, 0.2), mls in (4, 0)
        N = D == 1 ? 400 : 64
        eq, pg, ts, _ = setup(D = D, N = N, rf = rf, R = 5.5, method = (:RK4, :MUSCL, 4), mood = :U2, mls = mls)
        r = evolve!(eq, pg, ts; tmax = D == 1 ? 1.0 : 0.5)
        @printf("  D=%d rf=%.1f MLS_order=%d  |dVU|=%.2e  [%.3f, %.3f]\n", D, rf, mls, abs(r.dVU), r.umin, r.umax)
    end
end
if occursin('c', part)
    println("(c) time history, 1-D N=400 rf=0.2 R=2.5, MUSCL2 VK (shock forms at t = 1/π ≈ 0.32)")
    eq, pg, ts, _ = setup(D = 1, N = 400, rf = 0.2, R = 2.5, method = (:RK2, :MUSCL, 2), limiter = :VK)
    let m0 = Harness.mass(pg)[1], t = 0.0
        for tend in 0.1:0.1:3.0
            while t < tend - 1e-12
                dt = min(0.2 * HyperCloud.get_time_step(pg, eq, ts.divergence_interpolator), tend - t); ts(eq, pg, t, dt); t += dt
            end
            @printf("  t=%.1f  ΔΣVU/ΣVU = %+.3e\n", t, (Harness.mass(pg)[1] - m0) / m0)
        end
    end
end
if occursin('d', part)
    println("(d) best weights: V minimising Σ_pairs |V_i c_ij + V_j c_ji|² (Σ V = |Ω|) vs Shepard volumes, 1-D N=400 R=2.5")
    for rf in (0.05, 0.2)
        eq, pg, ts, _ = setup(D = 1, N = 400, rf = rf, R = 2.5, method = (:RK2, :MUSCL, 2), limiter = :VK)
        n = pg.meta.N; I, J, C = coefficients(pg, 1)
        rows = Int[]; cols = Int[]; vals = Float64[]; r = Ref(0)
        seen = Set{Tuple{Int, Int}}(); idx = Dict((I[k], J[k]) => k for k in eachindex(I))
        for k in eachindex(I)
            i, j = I[k], J[k]; (min(i, j), max(i, j)) in seen && continue; push!(seen, (min(i, j), max(i, j)))
            kk = idx[(j, i)]; r[] += 1
            push!(rows, r[]); push!(cols, i); push!(vals, C[k][1]); push!(rows, r[]); push!(cols, j); push!(vals, C[kk][1])
        end
        Amat = sparse(rows, cols, vals, r[], n)
        # min |A V|² s.t. Σ V = 1: KKT
        K = [Amat' * Amat ones(n); ones(n)' 0]; sol = K \ [zeros(n); 1.0]; Vb = sol[1:n] .* n .* (1 / n)
        Vs = pg.core.volumes[1:n] ./ sum(pg.core.volumes[1:n])
        @printf("  rf=%.2f  defect with Shepard V: %.2e   with best V: %.2e   (min/max best V × N: %.3f / %.3f)\n", rf,
                norm(Amat * Vs) / norm(abs.(Amat) * Vs), norm(Amat * Vb) / norm(abs.(Amat) * abs.(Vb)), minimum(Vb) * n, maximum(Vb) * n)
        # drift with both weights
        u0 = [pg.rhos[i][1] for i in 1:n]
        res = evolve!(eq, pg, ts; tmax = 1.0, V = Vb)
        @printf("           drift after the shock: Shepard V %.2e, best V %.2e\n", abs(res.dVU), abs(res.dX))
    end
end
