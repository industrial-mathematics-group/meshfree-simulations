# Violations per conserved component (ρ, ρu, ρv, E) and their local/global range, per record file of diag_u2.jl.
using Serialization, Statistics, Printf
for f in ARGS
    d = deserialize(f); r = filter(x -> !x.plateau && x.t > 0, d.recs)
    cnt = [count(x -> x.m == m, r) for m in 1:4]
    q(v) = isempty(v) ? NaN : median(v)
    println(@sprintf("%-22s n=%5d  by component ρ/ρu/ρv/E: %s   median R/G per comp: %s", basename(f), length(r), join(cnt, "/"),
                     join([@sprintf("%.1e", q([x.R / x.G for x in r if x.m == m])) for m in 1:4], "/")))
end
