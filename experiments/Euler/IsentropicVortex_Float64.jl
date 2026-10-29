# IsentropicVortex_Float64: RK4SSPMUSCL4(none,AR,MLS4) in Float64 on the CPU (local, 16 threads), N = 50 … 400, to
# check whether the growth of its L2 error in Float32 (IsentropicVortex_Float32, 141² → 400²) is a precision floor.
# Same parameters as `tools/run_sim.jl IsentropicVortex` with SIM_T=Float64, so `@sim IsentropicVortex_Float64` finds
# every run as cached. Check: julia --project=env_headless tools/check_config.jl IsentropicVortex_Float64
module _IsentropicVortex_Float64_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "IsentropicVortex.jl"))
end

function get_config()
    config = _IsentropicVortex_Float64_Base.get_config()
    config.shared[:T] = Float64
    empty!(config.varied)
    config.varied[:Grid_Ns__all] = [50, 71, 100, 141, 200, 283, 400]
    config.defaults = Symbol.(["RK4SSPMUSCL4(none,AR,MLS4)"])
    push!(config.sources, @__FILE__)
    return config
end
