# Conservation defect as a statistic, attached by the post-processing pass of `run_all_simulations` (the configs
# of this folder list this file in `sources` and name `conservation_defect!` as `postprocess`).
#
# A regular statistic (`calc_stat`) only sees the states of one time step, so it cannot form m(t) − m(0) or the
# boundary flux; the post-processing pass has the whole run and its parameters. Stored statistics:
#   :conservation_defect        time series  (m(t) − m(0) − B(t)) / |m(t)|        (per component)
#   :final_conservation_defect  scalar       the same at t_end
#   :boundary_flux              time series  B(t), the exact net inflow through the boundary
# with m(t) = h^DS Σ_i u_i, the discrete mass of the `:mass` statistic of the hyperbolic preset (h = equivalent
# particle spacing of the data; this is the quantity the scheme conserves, not Σ V_i u_i with Shepard volumes),
# and B(t) = 0 on periodic domains, B(t) = (f(u_L) − f(u_R)) · n t for a 1-D Riemann problem with outflow
# boundaries whose waves stay inside the domain until t_end (f the flux of the run's equation). Other setups have
# no exact boundary flux: the statistics are not added (the function warns once per run).
using HyperCloud, StaticArrays

# Exact net inflow per unit time, or `nothing` if the setup has none.
function defect_boundary_rate(params, ::Val{M}, ::Type{T}) where {M, T}
    get(params, :Grid_periodic, false) === true && return zero(SVector{M, T})
    if get(params, :IC_name, :none) === :riemann && params[:PDE_D] == 1
        ctx = Dict{Symbol, Any}(:T => T)
        eq = HyperCloud.build_equation(HyperCloud.extract_namespace(params, :PDE), ctx)
        uL = HyperCloud.State{M, T}(Tuple(T.(params[:IC_uL])))
        uR = HyperCloud.State{M, T}(Tuple(T.(params[:IC_uR])))
        n = T(first(get(params, :IC_n, (1.0,))))
        f(u) = SVector{M, T}(HyperCloud.flux(eq, u)[1])    # x-component of the flux
        return (f(uL) - f(uR)) * sign(n)
    end
    return nothing
end

function conservation_defect!(data::HyperCloud.AbstractSimData{D, DS, M, T}) where {D, DS, M, T}
    rate = defect_boundary_rate(data.params, Val(M), T)
    if rate === nothing
        @warn "conservation_defect!: no exact boundary flux for this setup, defect not computed" data.params[:IC_name]
        return false
    end
    measure = prod(data.domain.spacing[1:DS])
    mass = [sum(u) * measure for u in data.u]
    t = data.t
    flux = [rate * (t[n] - t[1]) for n in eachindex(t)]
    defect = [(mass[n] - mass[1] - flux[n]) ./ abs.(mass[n]) for n in eachindex(t)]
    add_stat!(data, :boundary_flux, flux, :time)
    add_stat!(data, :conservation_defect, defect, :time)
    add_stat!(data, :final_conservation_defect, last(defect), Symbol[])
    return true
end
