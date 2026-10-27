# clainshock with the old extrapolating (outflow) boundaries on all four sides: the control for the exact plane-shock
# boundary states of clainshock.jl (paper review 2026-10-26).
module _ClainShockOutflowBase
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "clainshock.jl"))
end

function get_config()
    config = _ClainShockOutflowBase.get_config()
    delete!(config.shared, :Grid_bc_map)
    config.shared[:Grid_bc] = _ClainShockOutflowBase.bc_map_outflow
    push!(config.sources, @__FILE__)
    return config
end
