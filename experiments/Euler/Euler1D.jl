function get_config()
    methods = get_master_method_dict()

    # Third order (MUSCL3 is third order since HyperCloudCore 0.6.2). At shocks a degree-3 divergence
    # (weights of both signs) is not robust: MOOD must lower the divergence degree with the order, so
    # these use the coupled MLS order (Scheme_MLS_order = 0); a fixed order 3 and the VK limiter blow up
    # on the 2-D problems (see Simulations REPORT_2026-10-18.md).
    for m in (create_method(:RK3, :MUSCL, 3; mood_criterion = :U2, mood_strategy = :EPD1, MOOD_delta_relax = 0.0,
                            Scheme_MLS_order = 0, ignore = HyperCloud.ignore_relax),
              create_method(:RK3, :MUSCL, 3; mood_criterion = :U2, mood_strategy = :EPD2, MOOD_delta_relax = 0.0,
                            Scheme_MLS_order = 0, ignore = HyperCloud.ignore_relax))
        push!(methods, m)
    end
    # fixed MLS order 3 (conservative mixing, no robust fallback) and VK, for comparison
    push!(methods, create_method(:RK3, :MUSCL, 3; mood_criterion = :U2, mood_strategy = :EPD1, MOOD_delta_relax = 0.0,
                              Scheme_MLS_order = 3, custom_label = "RK3MUSCL3(EPD1,MLS3)", ignore = HyperCloud.ignore_relax))
    push!(methods, create_method(:RK3, :MUSCL, 3; limiter = :VK, limiter_mode = :hard, Scheme_MLS_order = 3, ignore = HyperCloud.ignore_relax))
    # Assuming `sod_euler_params` and `ignore_relax` are defined as in your original script
    # push!(methods, ...) (Keep your method pushes here)
    
    varied_params = create_varied_dict()
    GAS_GAMMA_EULER = 1.4
    # Setup shared variables
    relax_vel, relax_idx = generate_relax_params(3.0, 1, 3)
    dummy_eq = EulerEquation(Val(1), Float64, GAS_GAMMA_EULER)
    
    # 1. Define Sod Shock Tube in Primitive Variables (rho, u, p)
    V_L = State{3}(1.0, 0.0, 1.0)
    V_R = State{3}(0.125, 0.0, 0.1)
    
    # 2. Convert to Conservative Variables (rho, m, E) for the solver
    U_L = prim2cons(dummy_eq, V_L)
    U_R = prim2cons(dummy_eq, V_R)
    
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 1,
         :PDE_gamma => GAS_GAMMA_EULER,
        
        :sim_func_name => :run_direct_simulation,
        :tmax => 0.2,
        :snapshots => 51,
        
        # Grid parameters mapping
        :Grid_domain => :rectangular,
        :Grid_mins => (-0.5,),
        :Grid_maxs => (0.5,),
        :Grid_Ns => (330,),
        :Grid_bc => Dict(1 => :fixed_dirichlet, 2 => :outflow),
        :Grid_randomness_factor => (0.2,),
        :Grid_seed => 10,
        :Grid_mover => :none,
        
        # Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 5.5,
        
        # Initial Condition mappings aligned to Riemann{D, M}
        :IC_name => :riemann,
        :IC_uL => Tuple(U_L),
        :IC_uR => Tuple(U_R),
        :IC_p0 => (0.0,),
        
        :Time_CFL => 0.2,
        :Scheme_MLS_order => 2,
        
        # Relaxation parameters
        :relax_velocities => relax_vel,
        :relax_indices => relax_idx, 
        :relax_epsilon => 1e-6
    )

    defaults = ["RK2Upwind1", "RK2MUSCL2(EPD1)", "RK2MUSCL2(VK)", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)", "RK3MUSCL3(EPD1,MLS3)", "RK3MUSCL3(VK)", "Analytical Solution"]
    
    return SimulationConfig(shared, methods, defaults; varied = varied_params, reference = "analytical_solution")
end
get_config()