isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))
# =========================================================================
# EXACT CONSTANTS FOR DOUBLE MACH REFLECTION
# =========================================================================
export euler_2d_DMR
# The conservative states: (rho, rho*u, rho*v, E)
# Pre-shock (Right): rho = 1.4, u = 0, v = 0, p = 1.0
const DMR_PRE_SHOCK  = State{4, Float64}(1.4, 0.0, 0.0, 2.5)

# Post-shock (Left): rho = 8.0, u = 7.1447, v = -4.125, p = 116.5
# Energy = 116.5 / 0.4 + 0.5 * 8.0 * (7.1447^2 + (-4.125)^2) = 563.5
const DMR_POST_SHOCK = State{4, Float64}(8.0, 57.15767664977295, -33.0, 563.5)


# =========================================================================
# CUSTOM BOUNDARY CONDITION FUNCTORS
# =========================================================================

# The conditions run as kernels (`foreach_particle`), so they work on every backend (CPU and GPU) and in
# the precision of the run: the states are converted to the element type of the buffer.

# 1. TOP WALL: Time-Dependent Exact Shock Condition
struct DMRTopBC <: AbstractBoundaryCondition end

function (::DMRTopBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    S = eltype(rhos_buffer)
    foreach_particle(_dmr_top!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.positions, pg.core.tags, pg.core.is_boundary, tag,
                      eltype(S)(t), S(DMR_POST_SHOCK), S(DMR_PRE_SHOCK)))
    return nothing
end

@inline function _dmr_top!(i, buf, pos, tags, is_boundary, tag, t, post, pre)
    if is_boundary[i] && tags[i] == tag
        x, y = pos[i][1], pos[i][2]
        x_shock = oftype(t, 1 / 6) + (y + 20 * t) / sqrt(oftype(t, 3))
        buf[i] = x < x_shock ? post : pre
    end
    return nothing
end

# 2. BOTTOM WALL: Slip Wall & Post-Shock Inflow
struct DMRBottomBC <: AbstractBoundaryCondition end

function (::DMRBottomBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    # Step 1: Run the standard OutflowBC to populate the ghost cells
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t)
    S = eltype(rhos_buffer)
    foreach_particle(_dmr_bottom!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.positions, pg.core.tags, pg.core.is_boundary, tag, S(DMR_POST_SHOCK)))
    return nothing
end

@inline function _dmr_bottom!(i, buf, pos, tags, is_boundary, tag, post)
    if is_boundary[i] && tags[i] == tag
        if pos[i][1] < 1 / 6
            buf[i] = post
        else
            # Slip wall: flip Y-momentum!
            U = buf[i]
            buf[i] = typeof(U)(U[1], U[2], -U[3], U[4])
        end
    end
    return nothing
end

# 3. LEFT WALL: Fixed Inflow
struct DMRLeftBC <: AbstractBoundaryCondition end

function (::DMRLeftBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    S = eltype(rhos_buffer)
    foreach_particle(_dmr_left!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.tags, pg.core.is_boundary, tag, S(DMR_POST_SHOCK)))
    return nothing
end

@inline function _dmr_left!(i, buf, tags, is_boundary, tag, post)
    (is_boundary[i] && tags[i] == tag) && (buf[i] = post)
    return nothing
end


# =========================================================================
# MAIN SIMULATION SETUP
# =========================================================================

function euler_2d_DMR()
    p0 = (1.0/6.0, 0.0)
    n  = (sqrt(3.0), -1.0) 
    
    bc_map = Dict{Int, AbstractBoundaryCondition}(
        1 => DMRLeftBC(),    # Left Wall
        2 => OutflowBC(),    # Right Wall 
        3 => DMRBottomBC(),  # Bottom Wall
        4 => DMRTopBC(),     # Top Wall
        5 => OutflowBC()     # Corners 
    )

    # 1. Flat Namespace Shared Configuration
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 2,
        :PDE_gamma => 1.4,

        :T => Float32,

        :sim_func_name => :run_direct_simulation,                
        :tmax => 0.2,
        :snapshots => 20,
        
        :Grid_domain => :rectangular,
        :Grid_mins => (0.0, 0.0),
        :Grid_maxs => (3.2, 1.0),
        :Grid_Ns => (1920, 600),
        :Grid_periodic => false,
        :Grid_bc_map => bc_map,
        :Grid_randomness_factor => (0.2, 0.2),
        :Grid_seed => 42,
        :Grid_mover => :none,
        
        # New Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,
        
        :IC_name => :riemann,
        :IC_uL => DMR_POST_SHOCK,
        :IC_uR => DMR_PRE_SHOCK,
        :IC_p0 => p0,
        :IC_n => n,


        :Time_CFL => 0.2,
        
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
    # MUSCL4 with the MLS (divergence) order coupled to the MOOD order (0) or fixed (4: conservative, but the
    # fallback is not first order; 2: conservative with the degree-1 divergence of MUSCL2)
    for (o, strategy, mls) in ((4, :EPD1, 0), (4, :EPD2, 0), (4, :EPD1, 4), (4, :EPD1, 2), (3, :EPD1, 2))
        push!(methods, create_method(Symbol(:RK, o), :MUSCL, o; mood_criterion = :U2, mood_strategy = strategy,
                                     MOOD_delta_relax = 0.0, Scheme_MLS_order = mls,
                                     custom_label = "RK$(o)MUSCL$(o)($strategy,MLS$mls)", ignore = HyperCloud.ignore_relax))
    end
    defaults = ["RK2MUSCL2(EPD1)", "RK2MUSCL2(EPD2)", "RK2MUSCL2(VK)", "RK2MUSCL1", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)"]
    return SimulationConfig(shared, methods, defaults; sources = [@__FILE__])
end
