function get_config()
    relax_vel, relax_idx = generate_relax_params(2.0, 1, 1)
    
    varied_params = create_varied_dict()
    push!(varied_params, :Grid_Ns__1 => [ceil.(Int,N) for N = 10 .^(1.5:.1:3.5)])
    
    shared = Dict{Symbol, Any}(
        :PDE_name => :linear,
        :PDE_velocities => ((1.0,),),
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
        :Grid_periodic => true,
        :Grid_randomness_factor => (0.0,),
        :Grid_seed => 42,
        :Grid_mover => :none,
        :Grid_remove_ghosts => false,
        
        # Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 5.5,
        
        # Initial Condition mappings aligned to Gauss{D, M}
        :IC_name => :gauss,
        :IC_a => (1.0,),          # Amplitude (Tuple wrapper for State{1})
        :IC_b => (0.0,),          # Center point
        :IC_width => 1.0,         # Pulse width
        
        :Time_CFL => 0.2,
        :Scheme_MLS_order => 5,
        
        # Relaxation parameters
        :relax_velocities => relax_vel, 
        :relax_indices => relax_idx,
        :relax_epsilon => 1e-6
    )

    base_methods = get_master_method_dict()
    defaults = ["RK4MUSCL2", "RK4MUSCL3", "RK4MUSCL4", "RK4MUSCL5"]
    
    return SimulationConfig(shared, base_methods, defaults; varied = varied_params, reference = "analytical_solution")
end