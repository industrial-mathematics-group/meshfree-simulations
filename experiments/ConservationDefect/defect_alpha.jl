isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# Weight: the parameter α of the exponential weight for MUSCL4(EPD1, k_div = 4).
function get_config()
    return defect_config(defect_shock_shared(), [(4, :EPD1, :AR, 4)]; varied = Dict(:Weight_alpha => [0.5, 1.0, 2.0, 4.0]))
end
