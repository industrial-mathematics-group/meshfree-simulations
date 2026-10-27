isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))
function get_config()
    relax_vel, relax_idx = generate_relax_params(2.0, 1, 1)
    
    varied_params = create_varied_dict()
    push!(varied_params, :Grid_Ns__all => [ceil.(Int,N) for N = 10 .^(1.5:.1:3.5)])
    
    shared = Dict{Symbol, Any}(
        :PDE_name => :burgers,
        :PDE_D => 1,
        :T => Float64,            # the sweeps of 2026-10-22 ran in Float64 (cluster, see reports/)
        
        :sim_func_name => :run_direct_simulation,
        :tmax => 10.0,
        :snapshots => 5,
        
        # Grid parameters mapping
        :Grid_domain => :rectangular,
        :Grid_mins => (-5.0,),
        :Grid_maxs => (5.0,),
        :Grid_Ns => (500,),
        :Grid_bc => Dict(0 => :outflow, 1 => :outflow),
        :Grid_periodic => false,
        :Grid_randomness_factor => (0.4,),
        :Grid_seed => 42,
        :Grid_mover => :none,
        :Grid_remove_ghosts => false,
        
        # Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 5.5,
        
        # Initial Condition Riemann
        :IC_name => :riemann,
        :IC_uL => (1.0,),       
        :IC_uR => (0.0,), 
        :IC_p0 => (-2.0,), 
        :IC_n0 => (1.,),
        
        :Time_CFL => 0.2,
        :Scheme_MLS_order => 4,
        
        # Relaxation parameters
        :relax_velocities => relax_vel, 
        :relax_indices => relax_idx,
        :relax_epsilon => 1e-6
    )

    base_methods = get_master_method_dict()

    add_review_methods!(base_methods)   # SSP steppers, AlignedRusanov flux, MLS 0/2/K (2026-10-26 setup)
    defaults = ["RK4MUSCL4(EPD1)", "RK4MUSCL4(EPD2)"]
    
    return SimulationConfig(shared, base_methods, defaults; varied = varied_params, reference = "analytical_solution")
end