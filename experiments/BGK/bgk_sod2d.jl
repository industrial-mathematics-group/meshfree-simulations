# Planar Sod shock tube for the 2-D discrete-velocity BGK model (study/bgk-lagrangian-contact): the 1-D Sod of bgk_sod.jl
# in the strip [-0.5, 0.5] × [0, 0.05], periodic in y, so the exact 1-D Euler solution applies in 2-D. τ = 1e-5, γ = 1.4
# (Chu-reduced pair, 8 × 8 Gauss–Hermite velocities as in bgk_vortex.jl), jitter 0.2, R = 3.5 (R = 2.5 failed on the
# stretched Lagrangian particles of the 1-D Sod), CFL 0.2, t = 0.2; Dirichlet left, outflow right. N_y = N_x / 20.
# Methods: ARS222MUSCL2(VK), ARS222MUSCL2(EPD2), ARS233MUSCL3(EPD2) of bgk_sod.jl, each static and Lagrangian (MLS 2 volumes,
# no particle management). The ALE flux test patch runs (tools/patches/bgk_ale_flux.jl) go to a separate save path.
#   SIM_BACKEND=cpu SIM_T=Float64 SIM_SAVE=data/sod2d SIM_SET="Grid_Ns=(200, 10)" julia --project=env_web -t 12 tools/run_sim.jl bgk_sod2d
module _bgk_sod2d_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "bgk_sod.jl"))
end

function get_config()
    shared = Dict{Symbol, Any}(
        :PDE_name => :bgk, :PDE_D => 2, :PDE_gamma => 1.4,
        :PDE_Nv => 8, :PDE_vmax => 7.0, :PDE_velocity_grid => :gauss_hermite,
        :PDE_tau => 1e-5,
        :sim_func_name => :run_bgk_simulation,
        :T => Float64, :tmax => 0.2, :snapshots => 10,
        :Grid_domain => :rectangular, :Grid_mins => (-0.5, 0.0), :Grid_maxs => (0.5, 0.05), :Grid_Ns => (200, 10),
        :Grid_periodic => (false, true),
        :Grid_bc_map => Dict{Int, AbstractBoundaryCondition}(1 => FixedDirichlet(), 2 => OutflowBC()),
        :Grid_randomness_factor => (0.2, 0.2), :Grid_seed => 10, :Grid_mover => :none,
        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 3.5,
        :IC_name => :riemann, :IC_uL => (1.0, 0.0, 0.0, 2.5), :IC_uR => (0.125, 0.0, 0.0, 0.25),
        :IC_p0 => (0.0, 0.0), :IC_n => (1.0, 0.0),
        :Time_CFL => 0.2,
    )
    base = _bgk_sod2d_Base.get_config().methods
    methods = Dict{Symbol, Dict{Symbol, Any}}()
    defaults = Symbol[]
    lag = Dict{Symbol, Any}(:Grid_mover => :lagrangian, :Grid_volumes => :mls, :Grid_volumes_degree => 2,
                            :Adapt_enabled => false)
    for name in ("ARS222MUSCL2(VK)", "ARS222MUSCL2(EPD2)", "ARS233MUSCL3(EPD2)")
        m = base[Symbol(name)]
        methods[Symbol("$name static")] = merge(m, Dict{Symbol, Any}(:Grid_mover => :none))
        methods[Symbol("$name lagrangian")] = merge(m, lag)
        push!(defaults, Symbol("$name static"), Symbol("$name lagrangian"))
    end
    return SimulationConfig(shared, methods, defaults; sources = [@__FILE__])
end
