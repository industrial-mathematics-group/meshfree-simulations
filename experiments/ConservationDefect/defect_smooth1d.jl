isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# Smooth periodic data in 1-D: advection of a Gaussian (experiments/Convergence/Convergence1D.jl), t = 10,
# θ = 0.2, 0.4, unstabilized K = 2 … 5 (k_div = K). The exact mass is constant.
function get_config()
    shared = defect_shock_shared()
    for k in (:Grid_bc, :IC_uL, :IC_uR, :IC_p0, :IC_n0); delete!(shared, k); end
    merge!(shared, Dict{Symbol, Any}(:PDE_name => :linear, :PDE_velocities => ((1.0,),), :Grid_periodic => true,
                                     :IC_name => :gauss, :IC_a => (1.0,), :IC_b => (0.0,), :IC_width => 1.0))
    specs = [(K, :none, :AR, K) for K in 2:5]
    return defect_config(shared, specs; varied = Dict(:Grid_randomness_factor__all => [0.2, 0.4]))
end
