# E2: library form (:gfd) vs conservative face forms (:sym, :symc) in the standalone solver.
# Mass conservation, smooth accuracy (pre-shock, exact solution) and post-shock accuracy (1-D: fine
# Godunov reference).
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf, StaticArrays, LinearAlgebra
D = parse(Int, get(ENV, "D", "1"))
u0 = D == 1 ? (x -> 1 + 0.5 * sinpi(2x[1])) : (x -> 1 + 0.5 * sinpi(2x[1]) * sinpi(2x[2]))
Ns = D == 1 ? (50, 100, 200, 400, 800) : (32, 64, 128)
tp = D == 1 ? 1.0 : 0.5                                     # post-shock time
if D == 1
    const REF = burgers_fv_1d(u0, 40000, tp)
    refval(x) = (h = REF[1][2] - REF[1][1]; REF[2][clamp(round(Int, x[1] / h + 0.5), 1, length(REF[2]))])
else                                                        # regular 256² run of the library form (conservative there)
    const NR = 256
    eqr, pgr, _, _ = setup(D = 2, N = NR, rf = 0.0, R = 2.5)
    Pr = Pairs(pgr; deg = 1)
    urr, _ = solve(Pr, nothing, eqr, pgr, [u0(pgr.core.positions[i]) for i in 1:pgr.meta.N], tp; order = 2, limiter = :bj)
    const REF = Dict(Tuple(round.(Int, pgr.core.positions[i] .* NR)) .% NR => urr[i] for i in 1:pgr.meta.N)
    refval(x) = REF[Tuple(mod.(round.(Int, x .* NR), NR))]
end
mkpath(joinpath(@__DIR__, "out")); out = open(joinpath(@__DIR__, "out", "e2_D$D.tsv"), "w")
println(out, "form\torder\tlimiter\tR\trf\tN\tsmooth_L1\tdVU_smooth\tpost_L1\tdVU_post\tdU_post\tumin\tumax")
for R in (2.5, 5.5), rf in (0.0, 0.2), N in Ns, (form, rs) in ((:gfd, :axis), (:symc, :normal), (:symc1, :normal)), (order, lim) in ((1, :none), (2, :none), (2, :bj))
    (rf == 0.0 && form != :gfd) && continue
    (D == 2 && R > 3 && N > 64 && form === :symc1) && continue     # the constraint solve gets too dense
    eq, pg, ts, _ = setup(D = D, N = N, rf = rf, R = R)
    n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]; V = collect(pg.core.volumes[1:n])
    P = Pairs(pg; deg = 1); A = face_coefficients(P, form)
    u_init = [u0(xi) for xi in x]
    m(u) = sum(V .* u)
    # smooth phase (t = 0.2 < 1/π)
    us, _ = solve(P, A, eq, pg, u_init, 0.2; order, limiter = lim, riemann = rs)
    ue = [burgers_exact(u0, xi, 0.2) for xi in x]
    eL1s = sum(V .* abs.(us .- ue))
    # post-shock
    up, k = solve(P, A, eq, pg, u_init, tp; order, limiter = lim, riemann = rs)
    eL1p = D == 1 || rf == 0.0 || true ? sum(V .* abs.(up .- [refval(xi) for xi in x])) : NaN
    println(out, join((form, order, lim, R, rf, N, eL1s, (m(us) - m(u_init)) / m(u_init), eL1p,
                       (m(up) - m(u_init)) / m(u_init), (sum(up) - sum(u_init)) / sum(u_init), minimum(up), maximum(up)), '\t'))
    flush(out)
end
close(out)
