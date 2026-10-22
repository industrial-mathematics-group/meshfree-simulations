# Mach3Step_1200x400: exactly the cluster runs of 2026-10-21/22 at 1200x400 (REPORT_2026-10-21.md), with the parameters
# `tools/run_sim.jl` used (SIM_SET "Grid_Ns=(1200,400); Time_CFL=…", SIM_T), so the fetched results in data/
# are found as cached runs instead of being rerun. Only the merged parameters are hashed (method entries
# override shared ones, labels are not part of the hash), so methods that ran at another CFL or precision are
# copies of their entry with that value as a per-method override.
isdefined(@__MODULE__, :euler_2d_mach3step) || Base.include(@__MODULE__, joinpath(@__DIR__, "Mach3Step.jl"))

function get_config()
    config = euler_2d_mach3step()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32, :Grid_Ns => (1200, 400), :Time_CFL => 0.8))
    config.defaults = ["RK2MUSCL2(EPD1)", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)"]
    config.methods[Symbol("RK3MUSCL3(EPD1,MLS2) CFL0.4")] = merge(config.methods[Symbol("RK3MUSCL3(EPD1,MLS2)")], Dict{Symbol, Any}(:Time_CFL => 0.4))
    push!(config.defaults, Symbol("RK3MUSCL3(EPD1,MLS2) CFL0.4"))
    return config
end
