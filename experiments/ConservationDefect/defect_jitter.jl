isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# Irregularity: the perturbation θ of the lattice for the unstabilized K = 2, 4 and MUSCL2/4(EPD1) (k_div = K).
function get_config()
    specs = [(2, :none, :AR, 2), (4, :none, :AR, 4), (2, :EPD1, :AR, 2), (4, :EPD1, :AR, 4)]
    return defect_config(defect_shock_shared(), specs; varied = Dict(:Grid_randomness_factor__all => [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]))
end
