# Canonical Shu–Osher forms of the explicit tableaus (Ketcheson 2008): which r make α_r, v_r ≥ 0, and what does
# checking Y_i instead of the Euler substeps miss for SSPRK(5,4)?
using HyperCloudCore, LinearAlgebra, Printf
function canonical(bt, r)
    s = num_stages(bt)
    K = zeros(s + 1, s + 1)
    K[1:s, 1:s] .= bt.a_t; K[s + 1, 1:s] .= bt.b_t
    α = r * K / (I + r * K)
    v = (I - α) * ones(s + 1)
    return α, v
end
ssp_ok(bt, r) = (α_v = canonical(bt, r); minimum(α_v[1]) >= -1e-12 && minimum(α_v[2]) >= -1e-12)
function ssp_coefficient(bt)
    ssp_ok(bt, 1e-8) || return 0.0
    lo, hi = 1e-8, 10.0
    for _ in 1:200
        mid = (lo + hi) / 2
        ssp_ok(bt, mid) ? (lo = mid) : (hi = mid)
    end
    return lo
end
for (name, bt) in (("Euler", RK1_Euler_Tableau(Float64)), ("RK2 SSP", RK2_SSP_Tableau(Float64)), ("RK2 Ralston", RK2_Ralston_Tableau(Float64)),
                   ("RK3 SSP", RK3_SSP_Tableau(Float64)), ("RK4 classical", RK4_Classical_Tableau(Float64)), ("SSPRK(5,4)", RK4_SSP_Tableau(Float64)))
    C = ssp_coefficient(bt)
    @printf "%-14s SSP coefficient C = %.6f   Euler substep ratio 1/C = %s\n" name C (C > 0 ? @sprintf("%.6f", 1 / C) : "—")
    if C > 0
        α, v = canonical(bt, C)
        s = num_stages(bt)
        for i in 2:(s + 1)
            terms = join([@sprintf("%.4f·W%d", α[i, j], j) for j in 1:(i - 1) if α[i, j] > 1e-12], " + ")
            @printf "    %s = %.4f·Uⁿ + %s\n" (i <= s ? "Y$i" : "Uⁿ⁺¹") v[i] terms
        end
    end
end
# SSPRK(5,4): the stage the current MOOD check sees at checkpoint i is Y_i; express the substep W_{i-1} it should
# check through Y_i, U^n and earlier substeps: W_{i-1} = (Y_i - v_i U^n - Σ_{j<i-1} α_ij W_j) / α_{i,i-1}.
bt = RK4_SSP_Tableau(Float64); C = ssp_coefficient(bt); α, v = canonical(bt, C); s = num_stages(bt)
println("\nSSPRK(5,4): substeps used by later stages (j → stages that use W_j):")
for j in 1:s
    users = [i for i in (j + 1):(s + 1) if α[i, j] > 1e-12]
    println("    W$j used in ", join([i <= s ? "Y$i" : "Uⁿ⁺¹" for i in users], ", "))
end
