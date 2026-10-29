# clainshock_800x800: the runs of the HyperCloudCore 0.7.8 rerun (reports/2026-10-29-rerun-078.md), GPU (CUDA), Float32,
# 800x800, CFL 0.8, with exactly the parameters the job script used (cluster/run_method.sh), so `@sim clainshock_800x800` finds every run as cached.
# Check: julia --project=env_headless tools/check_config.jl clainshock_800x800
module _clainshock_800x800_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "clainshock.jl"))
end

function get_config()
    config = _clainshock_800x800_Base.get_config()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32, :Grid_Ns => (800, 800), :Time_CFL => 0.8))
    config.defaults = Symbol.(["RK2MUSCL1",
                              "RK2MUSCL2(EPD1,AR,MLS2)",
                              "RK2MUSCL2(EPD2,AR,MLS2)",
                              "RK2MUSCL2(VK,AR,MLS2)",
                              "RK2MUSCL2(EPD2,RU,MLS2)",
                              "RK3MUSCL3(EPD1,AR,MLS2)",
                              "RK3MUSCL3(EPD1,AR,MLS3)",
                              "RK3MUSCL3(EPD1,AR,MLS0)",
                              "RK3MUSCL3(EPD2,AR,MLS2)",
                              "RK3MUSCL3(EPD2,AR,MLS3)",
                              "RK3MUSCL3(EPD2,AR,MLS0)",
                              "RK3MUSCL3(EPD2,RU,MLS2)",
                              "RK4SSPMUSCL4(EPD1,AR,MLS2)",
                              "RK4SSPMUSCL4(EPD1,AR,MLS4)",
                              "RK4SSPMUSCL4(EPD1,AR,MLS0)",
                              "RK4SSPMUSCL4(EPD2,AR,MLS2)",
                              "RK4SSPMUSCL4(EPD2,AR,MLS4)",
                              "RK4SSPMUSCL4(EPD2,AR,MLS0)",
                              "RK4SSPMUSCL4(EPD2,RU,MLS2)"])
    push!(config.sources, @__FILE__)
    return config
end
