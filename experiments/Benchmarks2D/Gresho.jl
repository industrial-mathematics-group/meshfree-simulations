# Gresho: the classic Gresho vortex (p0 = 5, peak Mach ≈ 0.38), 200², t = 3 (about 2.4 rotations of the core), see
# gresho_common.jl; errors against the exact (steady) solution.
isdefined(@__MODULE__, :GreshoIC) || Base.include(@__MODULE__, joinpath(@__DIR__, "gresho_common.jl"))
get_config() = gresho_config(5.0, 3.0, (200, 200), @__FILE__)
