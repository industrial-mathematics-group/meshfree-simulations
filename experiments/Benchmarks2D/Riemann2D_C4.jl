# Riemann2D_C4: Lax–Liu configuration 4 (four shocks), see riemann2d_common.jl.
isdefined(@__MODULE__, :LAX_LIU) || Base.include(@__MODULE__, joinpath(@__DIR__, "riemann2d_common.jl"))
get_config() = riemann2d_config(4, @__FILE__)
