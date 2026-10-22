# DMR_1920x600: exactly the cluster runs of 2026-10-21/22 at 1920x600 (REPORT_2026-10-21.md), with the parameters
# `tools/run_sim.jl` used (SIM_SET "Grid_Ns=(1920,600); Time_CFL=…", SIM_T), so the fetched results in data/
# are found as cached runs instead of being rerun. Only the merged parameters are hashed (method entries
# override shared ones, labels are not part of the hash), so methods that ran at another CFL or precision are
# copies of their entry with that value as a per-method override.
isdefined(@__MODULE__, :euler_2d_DMR) || Base.include(@__MODULE__, joinpath(@__DIR__, "DMR.jl"))

function get_config()
    config = euler_2d_DMR()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32, :Grid_Ns => (1920, 600), :Time_CFL => 0.8))
    config.defaults = ["RK2MUSCL2(EPD1)", "RK2MUSCL2(VK)", "RK2MUSCL2(EPD2)", "RK2MUSCL2(EPD0)", "RK2MUSCL2(SEPD0)", "RK2MUSCL2(EPD1,U1)", "RK2MUSCL2(BJ)", "RK2Upwind1", "RK4MUSCL4(EPD1,MLS0)", "RK4MUSCL4(EPD2,MLS0)", "RK4MUSCL4(EPD1,MLS2)"]
    push!(config.sources,@__FILE__)
    return config
end
