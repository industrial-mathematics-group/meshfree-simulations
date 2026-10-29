# DMR_960x300: the runs of the HyperCloudCore 0.7.8 rerun (reports/2026-10-29-rerun-078.md), GPU (CUDA), Float32,
# 960x300, CFL 0.8, with exactly the parameters the job script used (cluster/run_method.sh), so `@sim DMR_960x300` finds every run as cached.
# Check: julia --project=env_headless tools/check_config.jl DMR_960x300
module _DMR_960x300_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "DMR.jl"))
end

function get_config()
    config = _DMR_960x300_Base.euler_2d_DMR()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32, :Grid_Ns => (960, 300), :Time_CFL => 0.8))
    config.defaults = Symbol.(["RK4SSPMUSCL4(EPD1,AR,MLS0)",
                              "RK4SSPMUSCL4(EPD2,AR,MLS0)"])
    push!(config.sources, @__FILE__)
    return config
end
