# Relative L1 difference of the final density between CFL values of the same run (same grid, so the same
# particles): julia --project=env_gpu tools/cfl_compare.jl NAME "0.2,0.4,..."   (SIM_SET without Time_CFL)
using HyperCloud, StaticArrays, Printf
set_target_module!(Main)
set_save_path!(get(ENV, "SIM_SAVE", joinpath(@__DIR__, "..", "results")))
HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")
include(joinpath(@__DIR__, "config.jl"))
cfls = parse.(Float64, split(ARGS[2], ","))
for m in copy(config.defaults)
    HyperCloud.PDEStudioCore.is_reference_method(m) && continue
    base = assemble_params(config.shared, config.methods, m)
    rows = String[]; ref = nothing
    for c in cfls
        p = copy(base); p[:Time_CFL] = c
        does_sim_data_exist(p, Val(:raw)) || (push!(rows, @sprintf("CFL=%.2g: missing", c)); continue)
        d = load_sim_data(p, :raw)
        ρ = [s[1] for s in d.u[end]]
        ref === nothing && (ref = ρ)
        push!(rows, @sprintf("CFL=%.2g: %.2e (%d steps)", c, sum(abs, ρ .- ref) / sum(abs, ref), Int(d.stats[:k_step][1])))
    end
    println(rpad(string(m), 18), join(rows, " | "))
end
