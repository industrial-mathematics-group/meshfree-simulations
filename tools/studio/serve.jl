# Serves PDEStudio (WGLMakie) for a config of this repository on port 9384, so a browser (yours, or headless Chromium
# via browser.js) can use the studio on the stored runs in data/. For scripted sessions it also
#   - writes the studio state to state.json (status line, figure size, bounding boxes of all widgets), and
#   - evaluates Julia code dropped into cmd.jl (result or error in cmd.out), e.g. to inspect or set widgets.
#
#     cd /workspace/Simulations
#     julia -t 4 --project=env_web tools/studio/serve.jl MachJet80
using PDEStudio, WGLMakie, HyperCloud, StaticArrays
using WGLMakie.Bonito
const DIR = @__DIR__
const ROOT = abspath(joinpath(DIR, "..", ".."))

PDEStudioCore.set_target_module!(Main)
HyperCloud.CONFIG_ROOT_OVERWRITE[] = ROOT
set_save_path!(ROOT)
set_stat_preset!("hyperbolic")
WGLMakie.activate!()

name = get(ARGS, 1, "MachJet80")
config = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
config isa SimulationConfig || (config = Base.invokelatest(Main.get_config))
set_mode!(:lagrangian)
set_sim_config!(config)
fig = launch_plotter()

w() = PDEStudio.manager.widgets
js(x::AbstractString) = "\"" * escape_string(x) * "\""
js(x::Number) = isfinite(x) ? string(x) : "null"
js(x::AbstractVector) = "[" * join(js.(x), ",") * "]"
function bboxes()
    out = String[]
    for (k, v) in w()
        hasproperty(v, :layoutobservables) || continue
        r = v.layoutobservables.computedbbox[]
        push!(out, js(string(k)) * ": " * js([r.origin..., r.widths...]))
    end
    return "{" * join(out, ", ") * "}"
end
function write_state()
    st = haskey(w(), :status) ? replace(string(w()[:status][]), '\n' => ' ') : ""
    s = "{\"status\": " * js(st) * ", \"fig\": " * js(collect(fig.scene.viewport[].widths)) * ", \"widgets\": " * bboxes() * "}"
    write(joinpath(DIR, "state.json.tmp"), s); mv(joinpath(DIR, "state.json.tmp"), joinpath(DIR, "state.json"); force = true)
end

server = Bonito.Server(App(() -> fig), "0.0.0.0", 9384)
@info "Serving $name on $(Bonito.online_url(server, "/"))"
cmd_file, out_file = joinpath(DIR, "cmd.jl"), joinpath(DIR, "cmd.out")
rm(cmd_file; force = true)
while true
    try
        write_state()
        if isfile(cmd_file)
            code = read(cmd_file, String); rm(cmd_file)
            res = try
                "ok: " * repr(Base.invokelatest(include_string, Main, code))
            catch e
                "error: " * sprint(showerror, e, catch_backtrace())
            end
            write(out_file * ".tmp", res); mv(out_file * ".tmp", out_file; force = true)
        end
    catch e
        @error "studio loop" exception = (e, catch_backtrace())
    end
    sleep(0.2)
end
