# E1: mass drift of the library scheme on periodic Burgers vs jitter, radius, weight, resolution, order.
# Writes study/conservation/out/e1_D<D>.tsv.
include(joinpath(@__DIR__, "Harness.jl")); using .Harness, Printf
D = parse(Int, get(ENV, "D", "1"))
mkpath(joinpath(@__DIR__, "out"))
out = open(joinpath(@__DIR__, "out", "e1_D$D.tsv"), "w")
println(out, "method\tlimiter\tweight\tR\trf\tN\tdVU\tdU\tsteps\tumin\tumax")
tmax = D == 1 ? 1.0 : 0.5          # past the shock (t_s = 1/π ≈ 0.32)
Ns = D == 1 ? (50, 100, 200, 400, 800) : (32, 64, 128)
cases = [((:RK2, :MUSCL, 1), :none), ((:RK2, :MUSCL, 2), :VK), ((:RK2, :MUSCL, 2), :minmod)]
for (method, lim) in cases, w in (:exponential, :shifted_exp, :wendland), R in (2.5, 3.5, 5.5), rf in (0.0, 0.05, 0.1, 0.2, 0.3), N in Ns
    (rf == 0.0 && w != :exponential) && continue
    eq, pg, ts, _ = setup(D = D, N = N, rf = rf, R = R, weight = w, method = method, limiter = lim)
    r = evolve!(eq, pg, ts; tmax = tmax)
    line = join((join(string.(method)), lim, w, R, rf, N, r.dVU, r.dU, r.steps, r.umin, r.umax), '\t')
    println(out, line); flush(out)
end
close(out)
