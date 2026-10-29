# IsentropicVortex_Float32: the runs of the HyperCloudCore 0.7.8 rerun (reports/2026-10-29-rerun-078.md), GPU (CUDA), Float32
# (the CPU Float32 runs of the same jobs have the same parameters and are not kept separately),
# N = 50 … 400, CFL from the base config, with exactly the parameters the job script used (cluster/run_sweep.sh), so `@sim IsentropicVortex_Float32` finds every run as cached.
# Check: julia --project=env_headless tools/check_config.jl IsentropicVortex_Float32
module _IsentropicVortex_Float32_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "IsentropicVortex.jl"))
end

function get_config()
    config = _IsentropicVortex_Float32_Base.get_config()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32))
    empty!(config.varied)
    config.varied[:Grid_Ns__all] = [50, 71, 100, 141, 200, 283, 400]
    config.defaults = Symbol.(["RK4SSPMUSCL4(EPD2,AR,MLS2)",
                              "RK4SSPMUSCL4(EPD2,AR,MLS4)",
                              "RK4SSPMUSCL4(none,AR,MLS4)",
                              "RK3MUSCL3(EPD2,AR,MLS2)",
                              "RK3MUSCL3(none,AR,MLS3)",
                              "RK2MUSCL2(EPD2,AR,MLS2)"])
    push!(config.sources, @__FILE__)
    return config
end
