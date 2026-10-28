isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# Support: the weight range κ (in Δx) for MUSCL2(EPD1) and MUSCL4(EPD1, k_div = 4), aligned and geometric flux.
function get_config()
    specs = [(K, :EPD1, f, K) for K in (2, 4) for f in (:AR, :RU)]
    return defect_config(defect_shock_shared(), specs; varied = Dict(:Weight_range => [3.5, 4.5, 5.5, 7.5, 9.5]))
end
