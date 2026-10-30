# Shared setup of the 2-D Euler benchmarks (experiments/Benchmarks2D, /workspace/2D-Euler-Benchmarks.pdf): the method
# set, default parameters, slip walls and the conversion of primitive states. Included by every config of this folder.
#
# Defaults: GPU single precision (`:T => Float32`, the target setting; failures are rerun in Float64), jitter 0.2,
# exponential weights with R = 2.5 (the range of the shock problems), CFL 0.4, frozen MLS 2.
# Methods: the combinations known to work on HyperCloudCore 0.7.8 (labels of experiments/methods_review.jl):
# AlignedRusanov flux with a frozen MLS divergence degree (2 or K) and MOOD U2 (EPD1, EPD2) for K = 2 … 4, the VK
# limiter at K = 2, and MUSCL1. The cluster runs use `cluster/run_method.sh <config> <key> NX NY CFL cuda Float32`,
# which sets exactly :T, :Grid_Ns and :Time_CFL, so `@sim <config>` finds them when these equal the config's values.
# Configs include this file with `isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, …)`.
isdefined(@__MODULE__, :add_review_methods!) || Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))
using HyperCloud, StaticArrays

const BENCH_METHODS = ["RK2MUSCL1",
                       "RK2MUSCL2(VK,AR,MLS2)", "RK2MUSCL2(EPD1,AR,MLS2)", "RK2MUSCL2(EPD2,AR,MLS2)",
                       "RK3MUSCL3(EPD1,AR,MLS2)", "RK3MUSCL3(EPD2,AR,MLS2)", "RK3MUSCL3(EPD1,AR,MLS3)", "RK3MUSCL3(EPD2,AR,MLS3)",
                       "RK4SSPMUSCL4(EPD1,AR,MLS2)", "RK4SSPMUSCL4(EPD2,AR,MLS2)", "RK4SSPMUSCL4(EPD1,AR,MLS4)", "RK4SSPMUSCL4(EPD2,AR,MLS4)"]
# The matching keys of cluster/run_method.sh (same order)
const BENCH_KEYS = ["m1", "k2.vk.ar.mls2", "k2.epd1.ar.mls2", "k2.epd2.ar.mls2", "k3.epd1.ar.mls2", "k3.epd2.ar.mls2",
                    "k3.epd1.ar.mls3", "k3.epd2.ar.mls3", "k4.epd1.ar.mls2", "k4.epd2.ar.mls2", "k4.epd1.ar.mls4", "k4.epd2.ar.mls4"]

bench_methods() = add_review_methods!(get_master_method_dict())

"Conserved state (ρ, ρu, ρv, E) of the primitive state (ρ, u, v, p)."
bench_cons(γ, ρ, u, v, p) = (Float64(ρ), Float64(ρ * u), Float64(ρ * v), Float64(p / (γ - 1) + ρ * (u^2 + v^2) / 2))

"Shared parameters of a 2-D Euler benchmark; the pairs override or add entries (`:Grid_Ns => (400, 400)`)."
function bench_shared(overrides::Pair...)
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 2,
        :PDE_gamma => 1.4,
        :T => Float32,

        :sim_func_name => :run_direct_simulation,
        :snapshots => 20,

        :Grid_domain => :rectangular,
        :Grid_periodic => false,
        :Grid_randomness_factor => (0.2, 0.2),
        :Grid_seed => 42,
        :Grid_mover => :none,

        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,

        :Time_CFL => 0.4,
        :Scheme_MLS_order => 2,
    )
    for (k, v) in overrides
        shared[k] = v
    end
    return shared
end

# `sources`: the config file(s); they include this file behind `isdefined(@__MODULE__, :BENCH_METHODS) || …`
bench_config(shared; reference = nothing, postprocess = nothing, sources = String[]) =
    SimulationConfig(shared, bench_methods(), BENCH_METHODS; reference = reference, postprocess = postprocess,
                     sources = sources)

# =========================================================================
# Slip walls (rectangular domain tags: 1 left, 2 right, 3 bottom, 4 top): zero-gradient extrapolation, then the
# normal momentum is reflected. Kernels, so they run on every backend and in the precision of the run.
# =========================================================================
@inline function _bench_flip!(i, buf, tags, is_boundary, tag, d)
    if is_boundary[i] && tags[i] == tag
        U = buf[i]
        buf[i] = typeof(U)(U[1], d == 1 ? -U[2] : U[2], d == 2 ? -U[3] : U[3], U[4])
    end
    return nothing
end

"Slip wall normal to x (left/right boundary)."
struct BenchWallX <: AbstractBoundaryCondition end
"Slip wall normal to y (bottom/top boundary)."
struct BenchWallY <: AbstractBoundaryCondition end
for (W, d) in ((BenchWallX, 1), (BenchWallY, 2))
    @eval function (::$W)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
        OutflowBC()(pg, rhos_buffer, tag, ts, eq, t)
        foreach_particle(_bench_flip!, get_backend(), pg.meta.N, (rhos_buffer, pg.core.tags, pg.core.is_boundary, tag, $d))
        return nothing
    end
end
bench_walls() = Dict{Int, AbstractBoundaryCondition}(1 => BenchWallX(), 2 => BenchWallX(), 3 => BenchWallY(), 4 => BenchWallY())
bench_outflow() = Dict{Int, AbstractBoundaryCondition}(1 => OutflowBC(), 2 => OutflowBC(), 3 => OutflowBC(), 4 => OutflowBC())

"Pressure of a conserved state (ρ, ρu, ρv, E)."
@inline bench_pressure(U, γ) = (γ - 1) * (U[4] - (U[2]^2 + U[3]^2) / (2U[1]))
