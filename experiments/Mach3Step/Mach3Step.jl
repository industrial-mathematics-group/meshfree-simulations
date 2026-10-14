export euler_2d_mach3step

# State: (rho, rho*u, rho*v, E)
# Mach 3 Inflow: rho=1.4, u=3.0, p=1.0. 
# Energy = 1.0/0.4 + 0.5 * 1.4 * 9.0 = 8.8
const MACH3_INFLOW = State{4, Float64}(1.4, 4.2, 0.0, 8.8)

# =========================================================================
# 1. CUSTOM BOUNDARY CONDITIONS
# =========================================================================

struct HorizontalSlipWall <: AbstractBoundaryCondition end
function (::HorizontalSlipWall)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t) # Extrapolate from interior first
    @inbounds for i in 1:pg.meta.N
        if pg.core.tags[i] == tag && pg.core.is_boundary[i]
            U = rhos_buffer[i]
            rhos_buffer[i] = State{4, Float64}(U[1], U[2], -U[3], U[4]) # Flip Y-momentum
        end
    end
    return nothing
end

struct VerticalSlipWall <: AbstractBoundaryCondition end
function (::VerticalSlipWall)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t) # Extrapolate from interior first
    @inbounds for i in 1:pg.meta.N
        if pg.core.tags[i] == tag && pg.core.is_boundary[i]
            U = rhos_buffer[i]
            rhos_buffer[i] = State{4, Float64}(U[1], -U[2], U[3], U[4]) # Flip X-momentum
        end
    end
    return nothing
end

# =========================================================================
# 2. CUSTOM DOMAIN CONSTRUCTOR
# =========================================================================

function get_mach3step_domain( 
    mins::NTuple{2, Real}, 
    maxs::NTuple{2, Real},
    bc_map::Dict{Int, AbstractBoundaryCondition}
)
    T = Float64
    
    mins_f = Space{2, T}(mins...)
    maxs_f = Space{2, T}(maxs...)
    
    step_x = T(0.6)
    step_y = T(0.2)
    
    # Scale a small tolerance based on domain size for boundary tagging
    eps_tol = T(1e-5) * (maxs_f[1] - mins_f[1])

    is_interior_func = (pos) -> begin
        x, y = pos[1], pos[2]
        if x < mins_f[1] || x > maxs_f[1] || y < mins_f[2] || y > maxs_f[2]
            return false
        end
        # Cut out the step logic!
        if x >= step_x && y <= step_y
            return false
        end
        return true
    end

    actual_tag_func = (pos) -> begin
        if is_interior_func(pos)
            return 0
        end
        x, y = pos[1], pos[2]
        
        if x < mins_f[1] - eps_tol; return 1       # Inflow
        elseif x > maxs_f[1] + eps_tol; return 2   # Outflow
        elseif y > maxs_f[2] + eps_tol; return 4   # Top Wall
        elseif x >= step_x - eps_tol && y <= step_y + eps_tol
            if x < step_x + eps_tol; return 5      # Step Front
            else; return 6                         # Step Top
            end
        elseif y < mins_f[2] - eps_tol; return 3   # Bottom Wall
        else; return 2                             # Fallback
        end
    end

    # Added SVector{2, Bool}(false, false) to satisfy the new GeometricDomain signature
    return GeometricDomain(Val(:Mach3Step), mins_f, maxs_f, SVector{2, Bool}(false, false), is_interior_func, actual_tag_func, bc_map)
end

# =========================================================================
# 3. MAIN SIMULATION SETUP
# =========================================================================

function euler_2d_mach3step()
    mins = (0.0, 0.0)
    maxs = (3.0, 1.0)

    bc_map = Dict{Int, AbstractBoundaryCondition}(
        1 => FixedDirichlet(),
        2 => OutflowBC(),
        3 => HorizontalSlipWall(),
        4 => HorizontalSlipWall(),
        5 => VerticalSlipWall(),
        6 => HorizontalSlipWall()
    )

    # 1. Flat Namespace Shared Configuration
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 2,
        :PDE_gamma => 1.4,

        :sim_func_name => :run_direct_simulation,                
        :tmax => 4.0,
        :snapshots => 40,
        
        # Grid parameters mapping
        :Grid_domain => get_mach3step_domain(mins, maxs, bc_map),
        :Grid_mins => mins,
        :Grid_maxs => maxs,
        :Grid_Ns => (1800, 600),
        :Grid_periodic => false,
        :Grid_randomness_factor => (0.2, 0.2),
        :Grid_seed => 42,
        :Grid_mover => :none,
        
        # Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,
        
        # Initial Condition mappings directly aligned to Box{D, M}
        :IC_name => :box,
        :IC_u_bg => MACH3_INFLOW,
        :IC_u_box => MACH3_INFLOW,
        :IC_mins => mins,
        :IC_maxs => maxs,

        :Time_CFL => 0.2,
        :MLS_order => 2
    )

    # 2. Modernized Method Configurations
    methods = get_master_method_dict()
    defaults = ["RK2MUSCL2(EPD1)","RK2MUSCL2(EPD2)","RK2MUSCL2(VK)","RK2Upwind1"]
    
    return SimulationConfig(shared, methods, defaults)
end