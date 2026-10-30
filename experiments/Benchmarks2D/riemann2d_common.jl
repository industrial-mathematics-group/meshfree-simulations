# 2-D Riemann problems of Lax & Liu (1998), numbering and states as in Kurganov & Tadmor (2002): four constant states
# in the quadrants of [0, 1]² around (0.5, 0.5), γ = 1.4, outflow (zero-gradient) boundaries. Lax–Liu quadrants:
# 1 upper right, 2 upper left, 3 lower left, 4 lower right; `:q_riemann` orders them (BL, BR, TL, TR) = (3, 4, 2, 1).
# Configuration 3 (four shocks) is experiments/ClainShock. Defaults: 800², CFL 0.4, Float32 (bench_common.jl).
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

# configuration => ((ρ, u, v, p) of quadrants 1, 2, 3, 4), t_end
const LAX_LIU = Dict(
    4  => (((1.1, 0.0, 0.0, 1.1), (0.5065, 0.8939, 0.0, 0.35), (1.1, 0.8939, 0.8939, 1.1), (0.5065, 0.0, 0.8939, 0.35)), 0.25),
    6  => (((1.0, 0.75, -0.5, 1.0), (2.0, 0.75, 0.5, 1.0), (1.0, -0.75, 0.5, 1.0), (3.0, -0.75, -0.5, 1.0)), 0.3),
    12 => (((0.5313, 0.0, 0.0, 0.4), (1.0, 0.7276, 0.0, 1.0), (0.8, 0.0, 0.0, 1.0), (1.0, 0.0, 0.7276, 1.0)), 0.25),
    15 => (((1.0, 0.1, -0.3, 1.0), (0.5197, -0.6259, -0.3, 0.4), (0.8, 0.1, -0.3, 0.4), (0.5313, 0.1, 0.4276, 0.4)), 0.2),
    17 => (((1.0, 0.0, -0.4, 1.0), (2.0, 0.0, -0.3, 1.0), (1.0625, 0.0, 0.2145, 0.4), (0.5197, 0.0, -1.1259, 0.4)), 0.3),
)

function riemann2d_config(c::Int, file)
    q, tmax = LAX_LIU[c]
    γ = 1.4
    states = Tuple(bench_cons(γ, s...) for s in (q[3], q[4], q[2], q[1]))   # (BL, BR, TL, TR)
    shared = bench_shared(:tmax => tmax, :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => (800, 800),
                          :Grid_bc_map => bench_outflow(), :IC_name => :q_riemann, :IC_states => states, :IC_p0 => (0.5, 0.5))
    return bench_config(shared; sources = [joinpath(@__DIR__, "riemann2d_common.jl"), file])
end
