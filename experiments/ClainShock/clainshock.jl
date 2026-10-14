# Format: (rho, u, v, p)
const clain_riemann_init = (
    (0.138,  1.206, 1.206, 0.029), # 1: BL (Omega_ll)
    (0.5323, 0.0,   1.206, 0.3),   # 2: BR (Omega_lr)
    (0.5323, 1.206, 0.0,   0.3),   # 3: TL (Omega_ul)
    (1.5,    0.0,   0.0,   1.5)    # 4: TR (Omega_ur)
)

# Format: (rho, p, u, v) -> Reordered dynamically to (rho, u, v, p)
const implosionInit_raw = (
    (0.8,    1.0, 0.0, 0.0),
    (1.0,    1.0, 0.0, 0.7276),
    (1.0,    1.0, 0.7276, 0.0),
    (0.5313, 0.4, 0.0, 0.0)
)
const implosionNewOrder = Tuple((item[1], item[3], item[4], item[2]) for item in implosionInit_raw)
GAS_GAMMA_EULER = 1.4
function convertQuadrantInit(quadrant_data::NTuple{4, NTuple{4, T}}) where T <: AbstractFloat
    function primToCons(rho, u, v, p)
        mx = rho * u
        my = rho * v
        E = p / (GAS_GAMMA_EULER - 1.0) + 0.5 * rho * (u^2 + v^2)
        return (rho, mx, my, E)
    end
    return (
        primToCons(quadrant_data[1]...),
        primToCons(quadrant_data[2]...),
        primToCons(quadrant_data[3]...),
        primToCons(quadrant_data[4]...)
    )
end
bc_map = Dict(
    1 => :outflow,    # Left Wall
    2 =>  :outflow,    # Right Wall 
    3 =>  :outflow,  # Bottom Wall
    4 =>  :outflow,     # Top Wall
    5 =>  :outflow,    # Corners 
)
# =========================================================================
# MAIN SIMULATION SETUP
# =========================================================================

function get_config()
    relax_vel, relax_idx = generate_relax_params(4.0, 2, 4)
    
    # 1. Choose your test case
    # test_ic = convertQuadrantInit(implosionNewOrder)
    test_ic = convertQuadrantInit(clain_riemann_init)
    
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
         :PDE_gamma => GAS_GAMMA_EULER,
        :PDE_D => 2,
        
        :sim_func_name => :run_direct_simulation,
        :tmax => 0.3,
        :snapshots => 20,
        
        # Grid parameters mapping
        :Grid_domain => :rectangular,
        :Grid_mins => (0.0, 0.0),
        :Grid_maxs => (1.0, 1.0),
        :Grid_Ns => (800, 800),
        :Grid_bc => bc_map,
        :Grid_randomness_factor => (0.2, 0.2),
        :Grid_seed => 42,
        :Grid_mover => :none,
        
        # Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,
        
        # Initial Condition mappings aligned to QRiemann{D, M, N}
        :IC_name => :q_riemann,
        :IC_states => test_ic,
        :IC_p0 => (0.5, 0.5),
        
        :Time_CFL => 0.1,
        :MLS_order => 2,
        
        # Relaxation & Numerics
        :relax_velocities => relax_vel,
        :relax_indices => relax_idx,
        :relax_epsilon => 1e-6,
        :merge_factor => 0.0
    )

    method_dict = get_master_method_dict()
    defaults = ["RK2MUSCL2(EPD1)","RK2MUSCL2(EPD2)","RK2MUSCL2(VK)","RK2Upwind1"]
    
    return SimulationConfig(shared, method_dict, defaults)
end
