# clainshock_outflow_review_800x800: the cluster runs of the 2026-10-26 setup (reports/2026-10-27-review-rerun.md) at 800x800, CFL 0.8,
# Float64, with exactly the parameters `cluster/run_method.sh` used, so the fetched results (data/review/euler) are
# found as cached runs.
module _clainshock_outflow_review_800x800_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "ClainShock", "clainshock_outflow.jl"))
end

function get_config()
    config = _clainshock_outflow_review_800x800_Base.get_config()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float64, :Grid_Ns => (800, 800), :Time_CFL => 0.8))
    config.defaults = ["RK2MUSCL2(EPD2,AR,MLS2)",
                       "RK3MUSCL3(EPD2,AR,MLS2)",
                       "RK4SSPMUSCL4(EPD2,AR,MLS2)"]
    push!(config.sources, @__FILE__)
    return config
end
