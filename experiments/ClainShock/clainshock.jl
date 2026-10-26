isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))
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
bc_map_outflow = Dict(
    1 => :outflow,    # Left Wall
    2 =>  :outflow,    # Right Wall 
    3 =>  :outflow,  # Bottom Wall
    4 =>  :outflow,     # Top Wall
    5 =>  :outflow,    # Corners 
)

# Exact boundary states (paper review 2026-10-26): every boundary of the box is crossed by one of the four plane
# shocks between the quadrant states, whose speed follows from the Rankine–Hugoniot mass condition of the two
# adjacent states, s = (ρ_a u_a − ρ_b u_b) / (ρ_a − ρ_b) along the boundary direction. Until a shock reaches a
# corner (t > 1.1 here, tmax = 0.3) the state on each boundary is exactly the quadrant state on either side of
# the moving interface, so the ghosts can be prescribed instead of extrapolated (outflow), which removes the
# boundary from the list of possible error sources. Conserved states (rho, m_x, m_y, E), quadrant order
# (BL, BR, TL, TR) as in `QRiemann`.
"Shock speed along direction `d` between conserved states `a` and `b` (velocity component m_d / rho)."
_shock_speed(a, b, d) = (a[1 + d] - b[1 + d]) / (a[1] - b[1])

struct ExactQuadrantBC{S, T} <: AbstractBoundaryCondition
    states::S          # (BL, BR, TL, TR)
    p0::NTuple{2, T}   # centre
    mins::NTuple{2, T} # box
    maxs::NTuple{2, T}
    s_left::T          # speed of the BL/TL interface (along y, on x = xmin)
    s_right::T         # BR/TR interface (along y, on x = xmax)
    s_bottom::T        # BL/BR interface (along x, on y = ymin)
    s_top::T           # TL/TR interface (along x, on y = ymax)
end
function ExactQuadrantBC(states::NTuple{4, NTuple{4, T}}, p0::NTuple{2, T}, mins::NTuple{2, T}, maxs::NTuple{2, T}) where {T}
    BL, BR, TL, TR = states
    ExactQuadrantBC(states, p0, mins, maxs, _shock_speed(BL, TL, 2), _shock_speed(BR, TR, 2), _shock_speed(BL, BR, 1), _shock_speed(TL, TR, 1))
end
function (bc::ExactQuadrantBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    S = eltype(rhos_buffer); T = eltype(S)
    foreach_particle(_exact_quadrant!, get_backend(), pg.meta.N,
                     (rhos_buffer, pg.core.positions, pg.core.tags, pg.core.is_boundary, tag, T(t),
                      S(bc.states[1]), S(bc.states[2]), S(bc.states[3]), S(bc.states[4]),
                      T(bc.p0[1]), T(bc.p0[2]), T(bc.s_left), T(bc.s_right), T(bc.s_bottom), T(bc.s_top),
                      T(bc.mins[1]), T(bc.maxs[1]), T(bc.mins[2]), T(bc.maxs[2])))
    return nothing
end
@inline function _exact_quadrant!(i, buf, pos, tags, is_boundary, tag, t, BL, BR, TL, TR, x0, y0, sl, sr, sb, st,
                                  xmin, xmax, ymin, ymax)
    if is_boundary[i] && tags[i] == tag
        x, y = pos[i][1], pos[i][2]
        if x < xmin            # left boundary: BL below the BL/TL shock, TL above
            buf[i] = y < y0 + sl * t ? BL : TL
        elseif x > xmax        # right boundary: BR / TR
            buf[i] = y < y0 + sr * t ? BR : TR
        elseif y < ymin        # bottom boundary: BL / BR
            buf[i] = x < x0 + sb * t ? BL : BR
        else                   # top boundary: TL / TR
            buf[i] = x < x0 + st * t ? TL : TR
        end
    end
    return nothing
end
exact_bc(states) = let bc = ExactQuadrantBC(states, (0.5, 0.5), (0.0, 0.0), (1.0, 1.0)); Dict{Int, AbstractBoundaryCondition}(1 => bc, 2 => bc, 3 => bc, 4 => bc, 5 => bc); end
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
        :Grid_bc_map => exact_bc(test_ic),   # exact plane-shock boundary states; the outflow control: clainshock_outflow.jl
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
        :Scheme_MLS_order => 2,
        
        # Relaxation & Numerics
        :relax_velocities => relax_vel,
        :relax_indices => relax_idx,
        :relax_epsilon => 1e-6,
        :merge_factor => 0.0
    )

    method_dict = get_master_method_dict()

    add_review_methods!(method_dict)   # SSP steppers, AlignedRusanov flux, MLS 0/2/K (2026-10-26 setup)

    # Third order (MUSCL3 is third order since HyperCloudCore 0.6.2). At shocks a degree-3 divergence
    # (weights of both signs) is not robust: MOOD must lower the divergence degree with the order, so
    # these use the coupled MLS order (Scheme_MLS_order = 0); a fixed order 3 and the VK limiter blow up
    # on the 2-D problems (see Simulations REPORT_2026-10-18.md).
    for m in (create_method(:RK3, :MUSCL, 3; mood_criterion = :U2, mood_strategy = :EPD1, MOOD_delta_relax = 0.0,
                            Scheme_MLS_order = 0, ignore = HyperCloud.ignore_relax),
              create_method(:RK3, :MUSCL, 3; mood_criterion = :U2, mood_strategy = :EPD2, MOOD_delta_relax = 0.0,
                            Scheme_MLS_order = 0, ignore = HyperCloud.ignore_relax))
        push!(method_dict, m)
    end
    # MLS (divergence) order variants as in Mach3Step/DMR: 0 couples the divergence degree to the MOOD order,
    # a fixed order keeps it under order drops (conservative, but the fallback is not first order)
    for (o, strategy, mls) in ((3, :EPD1, 3), (3, :EPD1, 2), (4, :EPD1, 0), (4, :EPD1, 4), (4, :EPD2, 0), (4, :EPD1, 2))
        push!(method_dict, create_method(Symbol(:RK, o), :MUSCL, o; mood_criterion = :U2, mood_strategy = strategy,
                                         MOOD_delta_relax = 0.0, Scheme_MLS_order = mls,
                                         custom_label = "RK$(o)MUSCL$(o)($strategy,MLS$(mls))", ignore = HyperCloud.ignore_relax))
    end
    defaults = ["RK2MUSCL2(EPD1)","RK2MUSCL2(EPD2)","RK2MUSCL2(VK)","RK2MUSCL1", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)"]
    
    return SimulationConfig(shared, method_dict, defaults)
end
