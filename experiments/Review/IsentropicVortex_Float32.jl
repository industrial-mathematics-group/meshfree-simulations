# The isentropic-vortex runs in single precision (reports/2026-10-27-review-rerun.md §4): IsentropicVortex.jl with T = Float32.
module _VortexF32Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "Euler", "IsentropicVortex.jl"))
end

function get_config()
    config = _VortexF32Base.get_config()
    config.shared[:T] = Float32
    push!(config.sources, @__FILE__)
    return config
end
