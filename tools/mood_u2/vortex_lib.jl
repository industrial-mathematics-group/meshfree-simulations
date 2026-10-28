# Which part of the u2 criterion flags the smooth isentropic vortex? Runs the review's vortex with
# RK<K>MUSCL<K>(EPD2,AR,MLS2) and counts, for every DMP violation, why u2 rejected it.
#   julia -t 8 --project=env_mood vortex_mood.jl N K [tmax] [count=1]
using HyperCloud, StaticArrays, Printf
using Logging; global_logger(ConsoleLogger(stderr, Logging.Warn))
N, K = parse(Int, ARGS[1]), parse(Int, ARGS[2])
tmax = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 10.0
COUNT = length(ARGS) >= 4 ? ARGS[4] == "1" : true

const γ, β, C = 1.4, 5.0, (5.0, 5.0)
function vortex_state(x, y)
    dx, dy = x - C[1], y - C[2]; r2 = dx^2 + dy^2; ex = exp((1 - r2) / 2)
    u = 1 - β / (2π) * ex * dy; v = 1 + β / (2π) * ex * dx
    temp = 1 - (γ - 1) * β^2 / (8γ * π^2) * exp(1 - r2); ρ = temp^(1 / (γ - 1)); p = ρ^γ
    return (ρ, ρ * u, ρ * v, p / (γ - 1) + ρ * (u^2 + v^2) / 2)
end
struct VortexIC{T} <: HyperCloud.InitialCondition end
(::VortexIC{T})(pos) where {T} = State{4, T}(T.(vortex_state(Float64(pos[1]), Float64(pos[2]))))
HyperCloud.build_ic(::Val{:vortex}, conf::Dict, ctx::Dict) = VortexIC{ctx[:T]}()

const CNT = Dict{Symbol, Any}(); const VARIANT = "lib"

_, method = create_method(:RK3, :MUSCL, K; main_flux = :AlignedRusanov, mood_criterion = :U2, mood_strategy = :EPD2,
                          MOOD_delta_relax = parse(Float64, get(ENV, "DELTA_RELAX", "1e-4")), Scheme_MLS_order = 2)
K == 4 && (method[:Time_stepper] = :RK4SSP)
p = merge(Dict{Symbol, Any}(
    :PDE_name => :euler, :PDE_D => 2, :PDE_gamma => γ, :T => Float64, :sim_func_name => :run_direct_simulation,
    :tmax => tmax, :snapshots => 1, :Grid_domain => :rectangular, :Grid_mins => (0.0, 0.0), :Grid_maxs => (10.0, 10.0),
    :Grid_Ns => (N, N), :Grid_periodic => true, :Grid_randomness_factor => (0.2, 0.2), :Grid_seed => 42, :Grid_mover => :none,
    :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 2.5, :IC_name => :vortex, :Time_CFL => 0.8),
    method)
haskey(ENV, "VORTEX_NOMOOD") && (p[:MOOD_criterion] = :none)
haskey(ENV, "LIB_PLATEAU") && (p[:MOOD_plateau] = parse(Float64, ENV["LIB_PLATEAU"]))
haskey(ENV, "LIB_RATIO") && (p[:MOOD_curvature_ratio] = parse(Float64, ENV["LIB_RATIO"]))
enable_progress_bar!(false)
sim = run_direct_simulation(p)
sim === nothing && error("run failed")
x, u = sim.x[end], sim.u[end]
err = sum(abs(u[i][1] - vortex_state(x[i][1], x[i][2])[1]) for i in eachindex(u)) / length(u) * 100
c = Dict{Symbol, Int}()
@printf "[delta_relax=%s variant=%s] N=%d K=%d t=%.1f: L1(ρ)=%.3e min ρ=%.4f\n" get(ENV, "DELTA_RELAX", "1e-4") VARIANT N K tmax err minimum(v -> v[1], u)
if COUNT
    println("  criterion calls $(c[:calls]), DMP violations $(c[:dmp_fail]) (per component ρ, ρu, ρv, E: $([a[] for a in DMPCOMP]))")
    println("  u2 accepted $(c[:accept]); rejected: order<3 $(c[:low_order]) | sign test on a violating comp $(c[:sign_violating]), on a non-violating comp $(c[:sign_other]) | ratio test on a violating comp $(c[:ratio_violating]), on a non-violating comp $(c[:ratio_other])")
    println("  would be accepted if only the DMP-violating components were tested: $(c[:would_accept_violating_only])")
    println("  DMP violations by distance to (5,5), bins of 0.25 (meaningful for short tmax): ", [a[] for a in RHIST])
    println("  density violations, overshoot/δ in decades [<1e-2, 1e-2, 1e-1, 1, 10, ...]: ", [a[] for a in EXHIST])
    println("  density violations, local range/δ in decades [<1e-2, 1e-2, 1e-1, 1, 10, ...]: ", [a[] for a in RGHIST])
end
