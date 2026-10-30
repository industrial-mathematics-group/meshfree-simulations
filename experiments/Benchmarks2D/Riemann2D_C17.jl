# Riemann2D_C17: Lax–Liu configuration 17 (rarefaction, shock and two contacts (jet-like structure)), see riemann2d_common.jl.
isdefined(@__MODULE__, :LAX_LIU) || Base.include(@__MODULE__, joinpath(@__DIR__, "riemann2d_common.jl"))
get_config() = riemann2d_config(17, @__FILE__)
