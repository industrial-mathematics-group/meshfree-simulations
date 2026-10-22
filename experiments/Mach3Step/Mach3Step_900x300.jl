# Mach3Step_900x300: exactly the cluster runs of 2026-10-21/22 at 900x300 (REPORT_2026-10-21.md), with the parameters
# `tools/run_sim.jl` used (SIM_SET "Grid_Ns=(900,300); Time_CFL=…", SIM_T), so the fetched results in data/
# are found as cached runs instead of being rerun. Only the merged parameters are hashed (method entries
# override shared ones, labels are not part of the hash), so methods that ran at another CFL or precision are
# copies of their entry with that value as a per-method override.
isdefined(@__MODULE__, :euler_2d_mach3step) || Base.include(@__MODULE__, joinpath(@__DIR__, "Mach3Step.jl"))

function get_config()
    config = euler_2d_mach3step()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32, :Grid_Ns => (900, 300), :Time_CFL => 0.8))
    config.defaults = ["RK2MUSCL2(EPD1)", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)", "RK3MUSCL3(EPD1,MLS2)", "RK4MUSCL4(EPD1,MLS0)", "RK4MUSCL4(EPD2,MLS0)", "RK4MUSCL4(EPD1,MLS2)"]
    config.methods[Symbol("RK2MUSCL2(EPD1) Float64")] = merge(config.methods[Symbol("RK2MUSCL2(EPD1)")], Dict{Symbol, Any}(:T => Float64))
    push!(config.defaults, Symbol("RK2MUSCL2(EPD1) Float64"))
    config.methods[Symbol("RK4MUSCL4(EPD1,MLS0) Float64")] = merge(config.methods[Symbol("RK4MUSCL4(EPD1,MLS0)")], Dict{Symbol, Any}(:T => Float64))
    push!(config.defaults, Symbol("RK4MUSCL4(EPD1,MLS0) Float64"))
    return config
end
