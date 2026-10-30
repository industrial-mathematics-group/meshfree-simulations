# Riemann2D_C12: Lax–Liu configuration 12 (two shocks and two contacts), see riemann2d_common.jl.
isdefined(@__MODULE__, :LAX_LIU) || Base.include(@__MODULE__, joinpath(@__DIR__, "riemann2d_common.jl"))
get_config() = riemann2d_config(12, @__FILE__)
