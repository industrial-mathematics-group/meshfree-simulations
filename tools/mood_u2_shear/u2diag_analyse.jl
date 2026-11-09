# Quantiles of the U2 features per record file of diag_u2.jl (non-plateau violations).
#     julia tools/mood_u2_shear/u2diag_analyse.jl a.jls b.jls ...
# Distributions of the U2 features per case (records of diag_u2.jl), non-plateau violations only.
using Serialization, Statistics, Printf
q(v, p) = isempty(v) ? NaN : quantile(v, p)
fmt(v) = isempty(v) ? "      -" : @sprintf("%.2e/%.2e/%.2e/%.2e", q(v, 0.1), q(v, 0.5), q(v, 0.9), q(v, 0.99))
for f in ARGS
    d = deserialize(f)
    r = filter(x -> !x.plateau && x.t > 0, d.recs)
    isempty(r) && (println(basename(f), ": no records"); continue)
    τ = [d.K >= 3 ? x.τ2 : x.τ1 for x in r]
    println(@sprintf("%-18s K=%d  n=%6d   quantiles 10/50/90/99 %%", basename(f), d.K, length(r)))
    println("   v/R         ", fmt([x.v / x.R for x in r]))
    println("   v/τ1        ", fmt([x.v / max(x.τ1, 1e-300) for x in r]))
    d.K >= 3 && println("   v/τ2        ", fmt([x.v / max(x.τ2, 1e-300) for x in r]))
    println("   τ1/R        ", fmt([x.τ1 / x.R for x in r]))
    d.K >= 3 && println("   τ2/R        ", fmt([x.τ2 / x.R for x in r]))
    println("   mis/dU      ", fmt([x.mis / max(x.dU, 1e-300) for x in r]))
    println("   R/G         ", fmt([x.R / x.G for x in r]))
    if d.K >= 3
        println(@sprintf("   sign fail %.1f%%  ratio fail %.1f%%  negligible %.1f%%", 100mean(x -> !x.negl && !x.sign_ok, r),
                         100mean(x -> !x.negl && !x.ratio_ok, r), 100mean(x -> x.negl, r)))
    end
end
