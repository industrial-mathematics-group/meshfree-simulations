# Reconstructs a SimulationConfig from stored runs, for viewing them in an environment where the code that defined
# the original config can't be loaded (e.g. PDEStudio 0.6 / PDEStudioCore 0.5 next to HyperCloud 0.6, which still
# needs PDEStudioCore 0.4). The parameters stored with each run are already normalised, so a config built from them
# generates tasks with the same hashes: every task finds its file, nothing is recomputed (the solver is never called).
#
#   config_from_data(; filter = Dict(), sweep = Dict(:PDE_tau => :PDE_tau, :Grid_Ns => :Grid_Ns__1), label = f)
#     filter  parameter => value pairs a run must match (on the stored, normalised values)
#     sweep   stored key => varied key of the config (e.g. Grid_Ns => Grid_Ns__1 sweeps the first component)
#     label   params of a method -> method name (default: the differing keys and values)
# Parameters equal in all matching runs become shared, the swept keys become varied, and every distinct combination
# of the other differing keys becomes one method. `check_config(config)` counts the tasks without a file.
using PDEStudio
const _PSC = PDEStudio.PDEStudioCore

function config_from_data(; filter = Dict{Symbol, Any}(), sweep = Dict(:Grid_Ns => :Grid_Ns__1), label = nothing)
    ps = Dict{Symbol, Any}[]
    for f in _PSC._data_files()
        p = _PSC.JLD2.load(f, "params")
        all(_PSC.val2str(get(p, k, nothing)) == _PSC.val2str(v) for (k, v) in filter) && push!(ps, p)
    end
    isempty(ps) && error("no stored run matches $filter")
    str(p, k) = _PSC.val2str(get(p, k, nothing))
    allkeys = union(keys.(ps)...)
    differing = Set(k for k in allkeys if length(unique(str(p, k) for p in ps)) > 1)
    shared = Dict{Symbol, Any}(k => ps[1][k] for k in allkeys if !(k in differing))
    for k in keys(sweep)   # a varied key needs a default among the parameters
        k in differing && (shared[k] = ps[1][k])
    end
    method_keys = sort([k for k in differing if !haskey(sweep, k)], by = string)
    varied = create_varied_dict()
    for (k, vk) in sweep
        k in differing || continue
        vals = unique(p[k] for p in ps if haskey(p, k))
        vals = vk === :Grid_Ns__1 ? sort(unique(first.(vals))) : sort(vals)
        push!(varied, vk => vals)
    end
    methods = Dict{Symbol, Dict{Symbol, Any}}()
    for p in ps
        m = Dict{Symbol, Any}(k => p[k] for k in method_keys if haskey(p, k))
        name = Symbol(label === nothing ? join(("$k=$(str(p, k))" for k in method_keys if haskey(p, k)), " ") : label(m))
        methods[name] = m
    end
    names_ = sort(collect(keys(methods)), by = string)
    return SimulationConfig(shared, methods, names_; varied = varied)
end

function check_config(config)
    tasks = generate_tasks(config)
    missing_tasks = count(t -> !does_sim_data_exist(t, Val(:raw)), tasks)
    return (tasks = length(tasks), missing = missing_tasks)
end
