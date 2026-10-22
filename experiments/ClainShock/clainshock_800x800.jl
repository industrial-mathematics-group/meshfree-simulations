# clainshock_800x800: exactly the cluster runs of 2026-10-22 at 800x800, CFL 0.8, Float32
# (reports/2026-10-22-convergence-clainshock.md), with the parameters `tools/run_sim.jl` used, so the fetched
# results in data/ are found as cached runs instead of being rerun. clainshock.jl defines `get_config` itself,
# so it is loaded into a submodule.
module _ClainShockBase
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "clainshock.jl"))
end

function get_config()
    config = _ClainShockBase.get_config()
    merge!(config.shared, Dict{Symbol, Any}(:T => Float32, :Grid_Ns => (800, 800), :Time_CFL => 0.8))
    config.defaults = ["RK2MUSCL2(EPD1)", "RK3MUSCL3(EPD1)", "RK3MUSCL3(EPD2)", "RK3MUSCL3(EPD1,MLS2)", "RK3MUSCL3(EPD1,MLS3)",
                       "RK4MUSCL4(EPD1,MLS0)", "RK4MUSCL4(EPD2,MLS0)", "RK4MUSCL4(EPD1,MLS2)", "RK4MUSCL4(EPD1,MLS4)"]
    return config
end
