# L2 error of the isentropic vortex against its exact solution, attached by the post-processing pass of
# `run_all_simulations` (IsentropicVortex.jl lists this file in `sources` and names `vortex_error!` as `postprocess`,
# so `@sim IsentropicVortex…` adds it to new and to cached runs).
#
# The vortex is advected with the background velocity (1, 1) through the periodic box, so the exact solution at time
# t is the initial condition shifted by t (1, 1): u(x, t) = u₀((x − t (1, 1)) mod L). At t = 10 (one period, the
# `tmax` of the config) this is the initial condition itself. Stored statistics (per component ρ, ρu, ρv, E):
#   :vortex_l2error        time series  ‖u_h(t) − u(·, t)‖_L2 at every snapshot
#   :final_vortex_l2error  scalar       the same at t_end (= against the initial condition for t_end = 10)
# with ‖e‖_L2 = sqrt(h^DS Σ_i |e_i|²) over the particles, the measure of the `:l2error` statistic of the hyperbolic
# preset (h = equivalent particle spacing of the data). The exact state is evaluated in Float64 at the particle
# positions of each snapshot (also on moving grids) and the error is formed in Float64, then stored in the precision
# of the run (statistics must have the data's type).
using HyperCloud, StaticArrays

"Exact vortex state (ρ, ρu, ρv, E) at `x` and time `t` in the periodic box `[mins, maxs]`."
function vortex_exact(x, t, mins, maxs)
    L = maxs .- mins
    ξ = ntuple(d -> mins[d] + mod(Float64(x[d]) - t - mins[d], L[d]), 2)    # back along the characteristic (1, 1)
    return SVector{4, Float64}(vortex_state(ξ[1], ξ[2]))
end

function vortex_error!(data::HyperCloud.AbstractSimData{D, DS, M, T}) where {D, DS, M, T}
    haskey(data.stats, :final_vortex_l2error) && return false            # already computed: do not rewrite the file
    p = data.params
    mins, maxs = Float64.(p[:Grid_mins]), Float64.(p[:Grid_maxs])
    measure = prod(Float64.(data.domain.spacing[1:DS]))
    err = map(eachindex(data.t)) do n
        t = Float64(data.t[n])
        s = sum(eachindex(data.u[n])) do i
            abs2.(SVector{4, Float64}(data.u[n][i]) .- vortex_exact(data.x[n][i], t, mins, maxs))
        end
        SVector{M, T}(sqrt.(s .* measure))                                # stored in the precision of the data
    end
    add_stat!(data, :vortex_l2error, err, :time)
    add_stat!(data, :final_vortex_l2error, last(err), Symbol[])
    return true
end
