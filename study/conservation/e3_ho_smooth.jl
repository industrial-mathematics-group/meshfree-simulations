# E3: smooth convergence (RK4) of the GFD form vs conservative face forms exact to degree p, 1-D.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
D = 1; R = parse(Float64, get(ENV, "R", "5.5")); rf = 0.2
u0 = x -> 1 + 0.5 * sinpi(2x[1])
mkpath(joinpath(@__DIR__, "out")); out = open(joinpath(@__DIR__, "out", "e3_R$(R).tsv"), "w")
println(out, "pde\torder\tform\tN\tL1\tdVU")
for pde in (:advection, :burgers), o in (2, 3, 4), N in (50, 100, 200, 400, 800)
    eq, pg, ts, _ = setup(D = D, N = N, rf = rf, R = R, pde = pde)
    n = pg.meta.N; x = [pg.core.positions[i] for i in 1:n]; V = collect(pg.core.volumes[1:n])
    P = PairsHO(pg; q = o - 1)
    tmax = pde === :advection ? 1.0 : 0.2
    exact = pde === :advection ? [u0(xi .- tmax) for xi in x] : [burgers_exact(u0, xi, tmax) for xi in x]
    ui = [u0(xi) for xi in x]
    for (form, A) in (("gfd", nothing), ("face_p1", face_coefficients_p(P, 1)), ("face_p$(o-1)", o > 2 ? face_coefficients_p(P, o - 1) : nothing))
        (A === nothing && form != "gfd") && continue
        u, _ = solve_ho(P, A, eq, pg, ui, tmax; cfl = 0.2)
        println(out, join((pde, o, form, N, sum(V .* abs.(u .- exact)), (sum(V .* u) - sum(V .* ui)) / sum(V .* ui)), '\t')); flush(out)
    end
end
close(out)
