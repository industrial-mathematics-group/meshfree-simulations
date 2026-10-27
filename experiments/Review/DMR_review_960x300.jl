# DMR_review_960x300: the cluster runs of the 2026-10-26 setup (reports/2026-10-27-review-rerun.md) at 960x300, CFL 0.8,
# Float64, with exactly the parameters `cluster/run_method.sh` used, so the fetched results (data/review/euler) are
# found as cached runs.
# The coupled K = 4 runs (MLS 0) ran at 960x300 (they were expected to be the slowest).
module _DMR_review_960x300_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "DMR", "DMR.jl"))
end

function get_config()
    config = _DMR_review_960x300_Base.euler_2d_DMR()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float64, :Grid_Ns => (960, 300), :Time_CFL => 0.8))
    config.defaults = ["RK4SSPMUSCL4(EPD1,AR,MLS0)",
                       "RK4SSPMUSCL4(EPD2,AR,MLS0)"]
    push!(config.sources, @__FILE__)
    return config
end
