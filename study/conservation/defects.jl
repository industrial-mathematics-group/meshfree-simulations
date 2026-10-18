# Size of the two conservation defects of the MLS divergence (pair antisymmetry with Shepard volumes, row
# sums) versus jitter, weight function, support radius and fit degree.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
const D = parse(Int, get(ENV, "D", "1"))
const N = D == 1 ? 200 : 40
@printf("D=%d N=%d\n%-12s %4s %4s %3s | %9s %9s | %9s\n", D, N, "weight", "R", "rf", "deg", "asym(V)", "asym(1)", "rowsum")
for deg in (1, 3), R in (2.5, 3.5, 5.5), w in WEIGHTS, rf in (0.01, 0.05, 0.2)
    (deg == 3 && R < 3.5) && continue
    _, pg, _, _ = setup(D = D, N = N, rf = rf, R = R, weight = w, alpha = 1.0)
    I, J, C = coefficients(pg, deg)
    V = pg.core.volumes[1:pg.meta.N]
    a = residual_split(pg, I, J, C, V); b = residual_split(pg, I, J, C, ones(pg.meta.N))
    @printf("%-12s %4.1f %4.2f %3d | %9.2e %9.2e | %9.2e\n", w, R, rf, deg, a.asym, b.asym, a.rowsum)
end
