# Does a relaxed u2 still bound discontinuities? 1-D Sod (Euler) and 2-D box advection, MUSCL K with EPD2/U2,
# frozen MLS 2, AlignedRusanov, jitter 0.2, κ = 2.5 (the review setup).
#   U2_VARIANT=abc DELTA_ABS=1e-4 julia -t 8 --project=env_mood shock_mood.jl CASE K N
using HyperCloud, StaticArrays, Printf
using Logging; global_logger(ConsoleLogger(stderr, Logging.Warn))
case, K, N = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
const COUNT = true
include("u2_override.jl")

stepper = K == 4 ? :RK4SSP : (K == 3 ? :RK3 : :RK2)
D = case == "sod" ? 1 : 2
dx = case == "sod" ? 1.0 / N : 1.0 / N
dr = haskey(ENV, "DELTA_ABS") ? parse(Float64, ENV["DELTA_ABS"]) / dx^D : 1e-4
_, method = create_method(stepper, :MUSCL, K; main_flux = :AlignedRusanov, mood_criterion = :U2, mood_strategy = :EPD2,
                          MOOD_delta_relax = dr, Scheme_MLS_order = 2)
haskey(ENV, "NOMOOD") && (method[:MOOD_criterion] = :none)
common = Dict{Symbol, Any}(:T => Float64, :sim_func_name => :run_direct_simulation, :snapshots => 1,
    :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 2.5, :Grid_seed => 42, :Grid_mover => :none,
    :Grid_domain => :rectangular)
p = if case == "sod"
    merge(common, Dict{Symbol, Any}(:PDE_name => :euler, :PDE_D => 1, :PDE_gamma => 1.4, :tmax => 0.2,
        :Grid_mins => (-0.5,), :Grid_maxs => (0.5,), :Grid_Ns => (N,), :Grid_periodic => false,
        :Grid_bc => Dict{Int, AbstractBoundaryCondition}(1 => FixedDirichlet(), 2 => FixedDirichlet()),
        :Grid_randomness_factor => (0.2,), :IC_name => :riemann, :IC_uL => (1.0, 0.0, 2.5), :IC_uR => (0.125, 0.0, 0.25),
        :IC_p0 => (0.0,), :IC_n => (1.0,), :Time_CFL => 0.4), method)
else
    merge(common, Dict{Symbol, Any}(:PDE_name => :linear, :PDE_D => 2, :PDE_velocities => ((1.0,), (1.0,)), :tmax => 0.5,
        :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => (N, N), :Grid_periodic => true,
        :Grid_randomness_factor => (0.2, 0.2), :IC_name => :box, :IC_u_bg => (0.0,), :IC_u_box => (1.0,),
        :IC_mins => (0.25, 0.25), :IC_maxs => (0.6, 0.6), :Time_CFL => 0.4), method)
end
enable_progress_bar!(false)
exact = analytical_solution(p)
sim = run_direct_simulation(p)
sim === nothing && error("run failed")
x, u, t = sim.x[end], sim.u[end], sim.t[end]
err = sum(i -> abs(exact(SVector(Tuple(x[i])..., t))[1] - u[i][1]), eachindex(u)) / length(u)
lo, hi = extrema(v -> v[1], u)
blo, bhi = case == "sod" ? (0.125, 1.0) : (0.0, 1.0)
c = Dict(k => v[] for (k, v) in CNT)
@printf "%-5s K=%d N=%d variant=%-4s δ=%s: L1=%.3e  overshoot=%+.2e undershoot=%+.2e  (DMP violations %d, u2 accepted %d)\n" case K N VARIANT (haskey(ENV, "NOMOOD") ? "no MOOD" : get(ENV, "DELTA_ABS", "1e-4·Δx^D")) err hi - bhi blo - lo c[:dmp_fail] c[:accept]
