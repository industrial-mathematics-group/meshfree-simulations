# Core's pipeline test 3 (1-D box, MUSCL order K, U2/EPD1, δ = 1e-4 Δx) for plateau/ratio settings.
using HyperCloud, StaticArrays, Printf
using Logging; global_logger(ConsoleLogger(stderr, Logging.Warn))
K = parse(Int, get(ARGS, 1, "3"))
for (pl, r) in ((0.0, 0.5), (1e-3, 0.5), (1e-3, 0.0))
    p = Dict{Symbol, Any}(:sim_func_name => :run_direct_simulation, :snapshots => 1, :tmax => 2.0,
        :Grid_domain => :rectangular, :Grid_mins => (-1.0,), :Grid_maxs => (1.0,), :Grid_Ns => (100,),
        :Grid_periodic => true, :Grid_seed => 42, :Grid_randomness_factor => (0.0,),
        :Weight_name => :exponential, :Weight_range => 2.5, :Weight_alpha => 1.0,
        :PDE_name => :linear, :PDE_velocities => ((1.0,),),
        :IC_name => :box, :IC_u_bg => (0.0,), :IC_u_box => (1.0,), :IC_mins => (-0.5,), :IC_maxs => (0.5,),
        :Time_stepper => :RK3, :Time_CFL => 0.4, :Scheme_MLS_order => 0,
        :Scheme_name => :MUSCL, :Flux_name => :Rusanov, :Scheme_order => K, :Limiter_name => :none,
        :MOOD_criterion => :U2, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4,
        :MOOD_plateau => pl, :MOOD_curvature_ratio => r)
    enable_progress_bar!(false)
    sim = run_direct_simulation(p)
    lo, hi = extrema(v -> v[1], sim.u[end])
    @printf "K=%d plateau=%-6g ratio=%-4g: min %+.3e  max-1 %+.3e\n" K pl r lo hi - 1
end
