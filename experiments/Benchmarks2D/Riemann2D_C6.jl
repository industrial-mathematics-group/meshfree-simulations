# Riemann2D_C6: Lax–Liu configuration 6 (four contact discontinuities (vortex sheets)), see riemann2d_common.jl.
isdefined(@__MODULE__, :LAX_LIU) || Base.include(@__MODULE__, joinpath(@__DIR__, "riemann2d_common.jl"))
get_config() = riemann2d_config(6, @__FILE__)
