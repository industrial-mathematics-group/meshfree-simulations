# Shared by the scripts of this study: exact Sod solution, loading a stored run of a config, particle spacings.
using Logging, Printf
using HyperCloud, StaticArrays
const ROOT = abspath(joinpath(@__DIR__, "..", ".."))
set_target_module!(Main); HyperCloud.CONFIG_ROOT_OVERWRITE[] = ROOT

# exact Sod (γ = 1.4, left (ρ, u, p) = (1, 0, 1), right (0.125, 0, 0.1), diaphragm at x = 0): density at x, t
const γ = 1.4; const pst = 0.30313017805064707; const ust = 0.92745262004895057
const ρL, pL, ρR, pR = 1.0, 1.0, 0.125, 0.1
const cL = sqrt(γ * pL / ρL)
const ρsL = ρL * (pst / pL)^(1 / γ); const csL = sqrt(γ * pst / ρsL)
const ρsR = ρR * ((pst / pR + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * pst / pR + 1))
const S = ust * ρsR / (ρsR - ρR)              # shock speed
function ρexact(x, t)
    ξ = x / t
    ξ <= -cL && return ρL
    if ξ <= ust - csL
        c = 2 / (γ + 1) * (cL + (γ - 1) / 2 * (0 - ξ))
        return ρL * (c / cL)^(2 / (γ - 1))
    end
    ξ <= ust && return ρsL
    ξ <= S && return ρsR
    return ρR
end

const _CONFIGS = Dict{String, Any}()
function get_study_config(name)
    get!(_CONFIGS, name) do
        c = Base.include(Main, HyperCloud.resolve_config_path(Main, name))
        c isa SimulationConfig ? c : Base.invokelatest(Main.get_config)
    end
end

"Stored run of `method` of config `name` with N particles from the save path `save` (relative to the repository)."
function load_run(name, save, method, N)
    c = get_study_config(name)
    set_save_path!(isabspath(save) ? save : joinpath(ROOT, save))
    p = HyperCloud.PDEStudioCore.assemble_params(c.shared, c.methods, Symbol(method)); p[:Grid_Ns] = (N,)
    return load_sim_data(p, :raw)
end

"Sorted positions, densities, spacings and states of snapshot k (default: last) of a 1-D run."
function snapshot(d, k = length(d.t))
    x = [q[1] for q in d.x[k]]; pm = sortperm(x); x = x[pm]; u = d.u[k][pm]
    h = [i == 1 ? x[2] - x[1] : i == length(x) ? x[end] - x[end - 1] : (x[i + 1] - x[i - 1]) / 2 for i in eachindex(x)]
    return (; x, ρ = first.(u), h, u)
end

"Position where ρ crosses (ρ*_R + ρ_R)/2 (shock) in the last snapshot."
function shock_position(s)
    lvl = (ρsR + ρR) / 2
    i = findlast(k -> s.ρ[k] > lvl, 1:length(s.ρ) - 1)
    return s.x[i] + (s.x[i + 1] - s.x[i]) * (s.ρ[i] - lvl) / (s.ρ[i] - s.ρ[i + 1])
end

"Binned density profile of the last snapshot of a 2-D planar run: bin centres, mean ρ (NaN for empty bins)."
function binned_profile(d, nx; k = length(d.t))
    xs = range(-0.5, 0.5, length = nx + 1); xc = (xs[1:end-1] .+ xs[2:end]) ./ 2
    s = zeros(nx); n = zeros(Int, nx)
    for (q, u) in zip(d.x[k], d.u[k])
        b = clamp(1 + floor(Int, (q[1] + 0.5) * nx), 1, nx); s[b] += u[1]; n[b] += 1
    end
    return collect(xc), [n[b] > 0 ? s[b] / n[b] : NaN for b in 1:nx]
end
