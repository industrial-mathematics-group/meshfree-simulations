# Checks that `@sim NAME` would find every run of an experiment as cached: generates the tasks of each config
# (exactly as the pipeline does), looks each one up in the save path by its parameter hash, and lists the tasks
# without a file. With --orphans it also lists the data files that no task of the given configs matches, with the
# parameters in which the closest task differs (to see why a file is not found).
#
#     julia --project=env_headless tools/check_config.jl [--save DIR] [--orphans] NAME [NAME...]
#
# DIR defaults to the repository root (the save path of `startup.jl local`); the files are read from DIR/data.
using HyperCloud, StaticArrays, Printf
const PSC = HyperCloud.PDEStudioCore
const JLD2 = PSC.JLD2

args = copy(ARGS)
save = "."; orphans = false
while !isempty(args) && startswith(args[1], "--")
    a = popfirst!(args)
    a == "--save" ? (global save = popfirst!(args)) : a == "--orphans" ? (global orphans = true) : error("unknown option $a")
end
set_target_module!(Main)
set_save_path!(abspath(save))
HyperCloud.CONFIG_ROOT_OVERWRITE[] = joinpath(@__DIR__, "..")

function load_config(name)
    path = HyperCloud.resolve_config_path(Main, name)
    c = Base.include(Main, path)
    c isa Function && (c = Base.invokelatest(c))
    c isa SimulationConfig || (c = Base.invokelatest(Main.get_config))
    return c
end

taskhash = Dict{String, Tuple{String, Dict}}()
allok = true
for name in args
    config = load_config(name)
    tasks = Base.invokelatest(generate_tasks, config)
    missing_tasks = [t for t in tasks if !does_sim_data_exist(t, Val(:raw))]
    for t in tasks
        taskhash[PSC._param_hash(t)] = (name, t)
    end
    @printf("%-40s %4d tasks, %4d found, %4d missing\n", name, length(tasks), length(tasks) - length(missing_tasks), length(missing_tasks))
    for t in missing_tasks
        println("    missing: ", get(t, :custom_label, "?"), "  N = ", t[:Grid_Ns], "  T = ", t[:T])
    end
    global allok &= isempty(missing_tasks)
end

if orphans
    files = PSC._data_files()
    orph = [f for f in files if !haskey(taskhash, something(PSC._filename_hash(basename(f)), ""))]
    println("\n$(length(files)) data files, $(length(orph)) not matched by these configs")
    tasks = [t for (_, t) in values(taskhash)]
    for f in orph
        p = PSC._normalize_params(JLD2.load(f, "params"))
        isempty(tasks) && (println("  ", basename(f)); continue)
        # the task with the fewest differing keys
        diffs(t) = (tn = PSC._normalize_params(t);
                    [k for k in union(keys(p), keys(tn)) if PSC.val2str(get(p, k, nothing)) != PSC.val2str(get(tn, k, nothing))])
        best = argmin(t -> length(diffs(t)), tasks)
        d = diffs(best)
        println("  ", basename(f)[1:min(end, 40)], "…  ", get(p, :custom_label, "?"), "  N = ", get(p, :Grid_Ns, "?"),
                "  closest task differs in: ", join(("$k ($(get(p, k, "–")) vs $(get(PSC._normalize_params(best), k, "–")))" for k in d[1:min(end, 6)]), ", "))
    end
end
exit(allok ? 0 : 1)
