isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# Smooth periodic data in 2-D (Table 7): Burgers with a Gaussian before the shock forms
# (experiments/Convergence/ConvergenceBurgers.jl), t = 0.5, N = 32 … 200 per side, θ = 0, 0.2, 0.4, κ = 3.5, 5.5,
# MUSCL2/4 unstabilized, EPD1 and VK (k_div = K). The exact mass is constant.
function get_config()
    shared = Dict{Symbol, Any}(
        :PDE_name => :burgers, :PDE_D => 2, :T => Float64, :sim_func_name => :run_direct_simulation,
        :tmax => 0.5, :snapshots => 5,
        :Grid_domain => :rectangular, :Grid_mins => (-5.0, -5.0), :Grid_maxs => (5.0, 5.0), :Grid_Ns => (100, 100),
        :Grid_periodic => true, :Grid_randomness_factor => (0.0, 0.0), :Grid_seed => 42, :Grid_mover => :none,
        :Grid_remove_ghosts => false,
        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 5.5,
        :IC_name => :gauss, :IC_a => (1.0,), :IC_b => (0.0, 0.0), :IC_width => 1.0,
        :Time_CFL => 0.2)
    specs = [(K, v, :AR, K) for K in (2, 4) for v in (:none, :EPD1, :VK)]
    return defect_config(shared, specs; Ns = DEFECT_N_2D,
                         varied = Dict(:Grid_randomness_factor__all => [0.0, 0.2, 0.4], :Weight_range => [3.5, 5.5]))
end
