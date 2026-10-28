isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# The shock test at θ = 0.4, κ = 5.5: MOOD (EPD1, EPD2) with the divergence degree frozen at 1 (k_div = 2),
# coupled to the MOOD order (k_div = 0) or frozen at K, with the aligned (AR) and the geometric (RU) interface
# flux; the unstabilized schemes and the VK / minmod limiters (k_div = K) for comparison.
function get_config()
    specs = [[(K, v, f, m) for K in 2:5 for v in (:EPD1, :EPD2) for f in (:AR, :RU) for m in unique((2, 0, K))];
             [(K, v, :AR, K) for K in 2:5 for v in (:none, :VK, :minmod)]]
    return defect_config(defect_shock_shared(), specs)
end
