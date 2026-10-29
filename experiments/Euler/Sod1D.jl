# Sod1D: Sod shock tube in 1-D (HyperCloud's euler_1d_sod: [-0.5, 0.5], jitter 0.2, R = 5.5, CFL 0.2, t = 0.2,
# Dirichlet inflow / outflow), CPU, Float64, N = 100 … 1000. The error statistics of the hyperbolic preset
# (:l1error, :l2error, relative variants) are computed against the exact Riemann solution (`reference =
# :analytical_solution`, evaluated in Float64).
# Methods: the combinations known to work on HyperCloudCore 0.7.8, i.e. the AlignedRusanov flux with a fixed (frozen)
# MLS divergence degree and MOOD (EPD1, EPD2), MLS 2 and MLS K, K = 2 … 4, and the VK limiter at K = 2; MUSCL1 as the
# first-order baseline. Labels as in experiments/methods_review.jl.
#   julia --project=env_headless -t 16 tools/run_sim.jl Sod1D        (SIM_BACKEND=cpu SIM_T=Float64 SIM_SAVE=.)
isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))

const SOD1D_METHODS = ["RK2MUSCL1",
                       "RK2MUSCL2(VK,AR,MLS2)", "RK2MUSCL2(EPD1,AR,MLS2)", "RK2MUSCL2(EPD2,AR,MLS2)",
                       "RK3MUSCL3(EPD1,AR,MLS2)", "RK3MUSCL3(EPD2,AR,MLS2)", "RK3MUSCL3(EPD1,AR,MLS3)", "RK3MUSCL3(EPD2,AR,MLS3)",
                       "RK4SSPMUSCL4(EPD1,AR,MLS2)", "RK4SSPMUSCL4(EPD2,AR,MLS2)", "RK4SSPMUSCL4(EPD1,AR,MLS4)", "RK4SSPMUSCL4(EPD2,AR,MLS4)"]

function get_config()
    shared = copy(euler_1d_sod().shared)          # the setup of HyperCloud's example, methods from the master dict
    shared[:T] = Float64
    methods = get_master_method_dict()
    add_review_methods!(methods)
    varied = create_varied_dict()
    push!(varied, :Grid_Ns__1 => [ceil(Int, 10^x) for x in 2.0:0.1:3.0])
    return SimulationConfig(shared, methods, SOD1D_METHODS; varied = varied, reference = :analytical_solution,
                            sources = [@__FILE__])
end
