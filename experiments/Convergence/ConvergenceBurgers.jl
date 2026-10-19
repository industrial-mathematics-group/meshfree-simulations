function get_config()
    relax_vel, relax_idx = generate_relax_params(2.0, 2, 1)

    varied_params = create_varied_dict()
    push!(varied_params, :Grid_Ns => [(N, N) for N in ceil.(Int, 10 .^ (1.5:.1:3.0))])
    
    shared = Dict{Symbol, Any}(
        :PDE_name => :burgers,
        :PDE_D => 2,
        
        :sim_func_name => :run_direct_simulation,
        :tmax => 0.5,              # before the shock: width 1 breaks at t = e^(1/2) / 2 ≈ 0.82
        :snapshots => 5,
        
        # Grid parameters mapping
        :Grid_domain => :rectangular,
        :Grid_mins => (-5.0, -5.0),
        :Grid_maxs => (5.0, 5.0),
        :Grid_Ns => (100, 1000),
        :Grid_periodic => true,
        :Grid_randomness_factor => (0.0, 0.0),
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
        :IC_b => (0.0, 0.0),      # Center point
        :IC_width => 1.0,         # Pulse width (1.0: smooth across the periodic boundary to 1e-11)
        
        :Time_CFL => 0.2,
        :Scheme_MLS_order => 0,
    )

    base_methods = get_master_method_dict()
    defaults = ["RK4MUSCL2", "RK4MUSCL3", "RK4MUSCL4", "RK4MUSCL5"]
    
    return SimulationConfig(shared, base_methods, defaults; varied = varied_params, reference = "analytical_solution")
end