# Mach3Step_review_900x300: the cluster runs of the 2026-10-26 setup (reports/2026-10-27-review-rerun.md) at 900x300, CFL 0.8,
# Float64, with exactly the parameters `cluster/run_method.sh` used, so the fetched results (data/review/euler) are
# found as cached runs.
# Not included: RK4SSPMUSCL4(EPD2,AR,MLS4) (exceeded the 300 min limit, not run to the end).
module _Mach3Step_review_900x300_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "Mach3Step", "Mach3Step.jl"))
end

function get_config()
    config = _Mach3Step_review_900x300_Base.euler_2d_mach3step()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float64, :Grid_Ns => (900, 300), :Time_CFL => 0.8))
    config.defaults = ["RK2MUSCL1",
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
                       "RK4SSPMUSCL4(EPD2,AR,MLS0)",
                       "RK4SSPMUSCL4(EPD2,RU,MLS2)"]
    push!(config.sources, @__FILE__)
    return config
end
