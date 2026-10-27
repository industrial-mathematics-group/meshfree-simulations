isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))
export euler_2d_mach3step

# State: (rho, rho*u, rho*v, E)
# Mach 3 Inflow: rho=1.4, u=3.0, p=1.0. 
# Energy = 1.0/0.4 + 0.5 * 1.4 * 9.0 = 8.8
const MACH3_INFLOW = State{4, Float64}(1.4, 4.2, 0.0, 8.8)

# =========================================================================
# 1. CUSTOM BOUNDARY CONDITIONS
# =========================================================================

# Kernels (`foreach_particle`), so the walls work on every backend and in the precision of the run.
@inline function _flip_momentum!(i, buf, tags, is_boundary, tag, d)
    if is_boundary[i] && tags[i] == tag
        U = buf[i]
        buf[i] = typeof(U)(U[1], d == 1 ? -U[2] : U[2], d == 2 ? -U[3] : U[3], U[4])
    end
    return nothing
end

struct HorizontalSlipWall <: AbstractBoundaryCondition end
function (::HorizontalSlipWall)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t) # Extrapolate from interior first
    foreach_particle(_flip_momentum!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.tags, pg.core.is_boundary, tag, 2)) # Flip Y-momentum
    return nothing
end

struct VerticalSlipWall <: AbstractBoundaryCondition end
function (::VerticalSlipWall)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t) # Extrapolate from interior first
    foreach_particle(_flip_momentum!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.tags, pg.core.is_boundary, tag, 1)) # Flip X-momentum
    return nothing
end

# =========================================================================
# 2. CUSTOM DOMAIN CONSTRUCTOR
# =========================================================================

function get_mach3step_domain(
    ::Type{T},
    mins::NTuple{2, Real},
    maxs::NTuple{2, Real},
    bc_map::Dict{Int, AbstractBoundaryCondition}
) where {T}

    mins_f = Space{2, T}(mins...)
    maxs_f = Space{2, T}(maxs...)
    
    step_x = T(0.6)
    step_y = T(0.2)
    
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

    # Ghost tags, exactly complementary to `is_interior_func` (no tolerance bands: a ghost that lies on a boundary
    # line, e.g. x = 0 or y = 1 up to rounding, belonged to the "fallback" outflow tag before, ~60 ghosts at
    # 1800x600; paper review 2026-10-26).
    actual_tag_func = (pos) -> begin
        if is_interior_func(pos)
            return 0
        end
        x, y = pos[1], pos[2]
        if x < mins_f[1]; return 1                 # Inflow
        elseif x > maxs_f[1]; return 2             # Outflow
        elseif y > maxs_f[2]; return 4             # Top Wall
        elseif x >= step_x && y <= step_y
            # Inside the step block: the nearest face decides (front x = step_x vs top y = step_y)
            return x - step_x < step_y - y ? 5 : 6  # Step Front : Step Top
        else; return 3                             # Bottom Wall (y <= mins)
        end
    end

    # Added SVector{2, Bool}(false, false) to satisfy the new GeometricDomain signature
    return GeometricDomain(Val(:Mach3Step), mins_f, maxs_f, SVector{2, Bool}(false, false), is_interior_func, actual_tag_func, bc_map)
end

# Area of the channel minus the step (needed for the volume normalisation).
HyperCloud.HyperCloudCore.domain_volume(geom::GeometricDomain{Val{:Mach3Step}, 2, T}) where {T} =
    prod(geom.maxs .- geom.mins) - (geom.maxs[1] - T(0.6)) * T(0.2)

# `:Grid_domain => :mach3step` builds the domain in the precision of the run (`:T`).
HyperCloud.build_domain(::Val{:mach3step}, conf::Dict, ctx::Dict) =
    get_mach3step_domain(real_type(ctx), conf[:mins], conf[:maxs], parse_bc(conf))

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
        :Grid_domain => :mach3step,
        :Grid_bc_map => bc_map,
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
        :Scheme_MLS_order => 2
    )

    # 2. Modernized Method Configurations
    methods = get_master_method_dict()
    add_review_methods!(methods)   # SSP steppers, AlignedRusanov flux, MLS 0/2/K (2026-10-26 setup)

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
    # MLS (divergence) order variants: 0 couples the divergence degree to the MOOD order (can lose mass),
    # a fixed order keeps it under order drops (conservative, but the fallback is no longer first order,
    # so it can reach negative pressure; PAD in MOOD rejects those candidates).
    for (o, strategy, mls) in ((3, :EPD1, 3), (3, :EPD1, 2), (4, :EPD1, 0), (4, :EPD1, 4), (4, :EPD2, 0), (4, :EPD1, 2))
        push!(methods, create_method(Symbol(:RK, o), :MUSCL, o; mood_criterion = :U2, mood_strategy = strategy,
                                     MOOD_delta_relax = 0.0, Scheme_MLS_order = mls,
                                     custom_label = "RK$(o)MUSCL$(o)($strategy,MLS$(mls))", ignore = HyperCloud.ignore_relax))
    end
    defaults = ["RK2MUSCL2(EPD1)","RK2MUSCL2(EPD2)","RK2MUSCL2(VK)","RK2MUSCL1", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)"]
    
    return SimulationConfig(shared, methods, defaults; sources = [@__FILE__])
end