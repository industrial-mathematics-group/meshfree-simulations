# MOOD substep checks (SSP) and the K = 1 orientation: one run per call.
#   julia -t 8 --project=env_ssp ssp_study.jl CASE K N SUBSTEPS [PLATEAU]
#   CASE: vortex | sod | box | r123 ; K: 1 (RK2, MUSCL1, no MOOD), 2 (RK2), 3 (RK3), 4 (SSPRK54); SUBSTEPS: 1 | 0
using HyperCloud, StaticArrays, Printf
using Logging; global_logger(ConsoleLogger(stderr, Logging.Warn))
case, K, N, sub = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]), ARGS[4] == "1"
plateau = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : nothing

const γ = 1.4
function vortex_state(x, y)
    β, c = 5.0, (5.0, 5.0)
    dx, dy = x - c[1], y - c[2]; r2 = dx^2 + dy^2; ex = exp((1 - r2) / 2)
    u = 1 - β / (2π) * ex * dy; v = 1 + β / (2π) * ex * dx
    temp = 1 - (γ - 1) * β^2 / (8γ * π^2) * exp(1 - r2); ρ = temp^(1 / (γ - 1)); p = ρ^γ
    return (ρ, ρ * u, ρ * v, p / (γ - 1) + ρ * (u^2 + v^2) / 2)
end
struct VortexIC{T} <: HyperCloud.InitialCondition end
(::VortexIC{T})(pos) where {T} = State{4, T}(T.(vortex_state(Float64(pos[1]), Float64(pos[2]))))
HyperCloud.build_ic(::Val{:vortex}, conf::Dict, ctx::Dict) = VortexIC{ctx[:T]}()

stepper = K == 4 ? :RK4SSP : (K == 3 ? :RK3 : :RK2)
_, method = create_method(stepper, :MUSCL, K; main_flux = :AlignedRusanov, mood_criterion = K == 1 ? :none : :U2,
                          mood_strategy = :EPD2, MOOD_delta_relax = 1e-4, Scheme_MLS_order = 2)
method[:Time_mood_substeps] = sub
plateau === nothing || (method[:MOOD_plateau] = plateau)
common = Dict{Symbol, Any}(:T => Float64, :sim_func_name => :run_direct_simulation, :snapshots => 1, :Grid_domain => :rectangular,
    :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 2.5, :Grid_seed => 42, :Grid_mover => :none)
dirichlet = Dict{Int, AbstractBoundaryCondition}(1 => FixedDirichlet(), 2 => FixedDirichlet())
p = if case == "vortex"
    merge(common, Dict{Symbol, Any}(:PDE_name => :euler, :PDE_D => 2, :PDE_gamma => γ, :tmax => 10.0,
        :Grid_mins => (0.0, 0.0), :Grid_maxs => (10.0, 10.0), :Grid_Ns => (N, N), :Grid_periodic => true,
        :Grid_randomness_factor => (0.2, 0.2), :IC_name => :vortex, :Time_CFL => 0.8))
elseif case in ("sod", "r123", "vac")
    uL, uR = case == "sod" ? ((1.0, 0.0, 2.5), (0.125, 0.0, 0.25)) : case == "r123" ? ((1.0, -2.0, 3.0), (1.0, 2.0, 3.0)) : ((1.0, -3.5, 7.125), (1.0, 3.5, 7.125))
    merge(common, Dict{Symbol, Any}(:PDE_name => :euler, :PDE_D => 1, :PDE_gamma => γ, :tmax => case == "sod" ? 0.2 : (case == "vac" ? 0.1 : 0.15),
        :Grid_mins => (-0.5,), :Grid_maxs => (0.5,), :Grid_Ns => (N,), :Grid_periodic => false, :Grid_bc => dirichlet,
        :Grid_randomness_factor => (0.2,), :IC_name => :riemann, :IC_uL => uL, :IC_uR => uR, :IC_p0 => (0.0,), :IC_n => (1.0,),
        :Time_CFL => parse(Float64, get(ENV, "CFL", "0.4"))))
else
    merge(common, Dict{Symbol, Any}(:PDE_name => :linear, :PDE_D => 2, :PDE_velocities => ((1.0,), (1.0,)), :tmax => 0.5,
        :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => (N, N), :Grid_periodic => true,
        :Grid_randomness_factor => (0.2, 0.2), :IC_name => :box, :IC_u_bg => (0.0,), :IC_u_box => (1.0,),
        :IC_mins => (0.25, 0.25), :IC_maxs => (0.6, 0.6), :Time_CFL => 0.4))
end
p = merge(p, method)
enable_progress_bar!(false)
t0 = time()
sim = run_direct_simulation(p)
wall = time() - t0
tag = @sprintf "%-6s K=%d N=%-4d substeps=%d%s" case K N sub (plateau === nothing ? "" : " plateau=$plateau")
if sim === nothing
    println(tag, ": FAILED (see log)")
    exit(0)
end
x, u, t = sim.x[end], sim.u[end], sim.t[end]
finite = all(v -> all(isfinite, v), u)
if case == "vortex"
    err = sum(i -> abs(u[i][1] - vortex_state(x[i][1], x[i][2])[1]), eachindex(u)) / length(u) * 100
    @printf "%s: L1(ρ)=%.4e min ρ=%.4f wall=%.1fs\n" tag err minimum(v -> v[1], u) wall
else
    ex = analytical_solution(p)
    ex === nothing && (ex = st -> SVector(NaN))
    err = sum(i -> abs(ex(SVector(Tuple(x[i])..., t))[1] - u[i][1]), eachindex(u)) / length(u)
    if case == "box"
        lo, hi = extrema(v -> v[1], u)
        @printf "%s: L1=%.4e overshoot=%+.2e undershoot=%+.2e wall=%.1fs\n" tag err hi - 1 -lo wall
    else
        ρmin = minimum(v -> v[1], u)
        pmin = minimum(v -> (γ - 1) * (v[3] - 0.5 * v[2]^2 / v[1]), u)
        ρmax = maximum(v -> v[1], u)
        @printf "%s: L1(ρ)=%.4e min ρ=%.3e min p=%.3e max ρ-1=%+.2e finite=%s wall=%.1fs\n" tag err ρmin pmin ρmax - 1 finite wall
    end
end
