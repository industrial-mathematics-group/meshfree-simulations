# Shared u2 override with counters and variant switches (U2_VARIANT letters a/b/c), see vortex_mood.jl.
const HCC = HyperCloud.HyperCloudCore
# --- counters (atomic: the criterion runs inside a threaded kernel) ---
const CNT = Dict(k => Threads.Atomic{Int}(0) for k in
    (:calls, :dmp_fail, :accept, :low_order, :sign_violating, :sign_other, :ratio_violating, :ratio_other, :would_accept_violating_only))
const DMPCOMP = [Threads.Atomic{Int}(0) for _ in 1:4]
const RHIST = [Threads.Atomic{Int}(0) for _ in 1:24]          # DMP violations by distance to (5, 5), bins of 0.25
const EXHIST = [Threads.Atomic{Int}(0) for _ in 1:12]         # density overshoot / δ, log10 bins: 1 = <1e-2, ..., 12 = ≥1e9
const RGHIST = [Threads.Atomic{Int}(0) for _ in 1:12]         # local density range / δ, same bins
inc(k) = Threads.atomic_add!(CNT[k], 1)
const VARIANT = get(ENV, "U2_VARIANT", "")              # letters: a = units, b = violating components only, c = no ratio test
const VARIANT_UNITS = occursin('a', VARIANT)
const VARIANT_COMPS = occursin('b', VARIANT)
const VARIANT_NORATIO = occursin('c', VARIANT)
const VARIANT_PLATEAU = occursin('p', VARIANT)
const ETA = parse(Float64, get(ENV, "U2_ETA", "1e-3"))
const GLO = Ref{Any}(nothing); const GHI = Ref{Any}(nothing)   # global extrema of U_prev per component, per MOOD check

# detect_and_drop! of main plus the global extrema of the previous stage (one reduction per check)
@eval HCC function detect_and_drop!(
    mood::MOOD, work::MOODWork, interp::DivergenceInterpolator, pg::ParticleGrid{D, M, T},
    rho_prev::AbstractVector, rho_cand::AbstractVector, full::Bool; eq = nothing
) where {D, M, T}
    mood.criterion isa NoMOOD && return false
    N = pg.meta.N
    Main.GLO[] = mapreduce(identity, (a, b) -> min.(a, b), view(rho_prev, 1:N))
    Main.GHI[] = mapreduce(identity, (a, b) -> max.(a, b), view(rho_prev, 1:N))
    nb_slices = pg.neighbor.ranges; is_boundary = pg.core.is_boundary; orders = pg.core.particle_orders
    triggered, changed, dirty, snapshot = parent(work.triggered), parent(work.changed), parent(work.dirty), parent(work.snapshot)
    triggered_step = parent(work.triggered_step)
    fill!(view(triggered, 1:N), false); fill!(view(changed, 1:N), false)
    foreach_particle(_mood_flag!, get_backend(), N, (mood, interp, pg, rho_prev, rho_cand, nb_slices, is_boundary, dirty, triggered, full, eq))
    foreach_particle(_mood_drop!, get_backend(), N, (orders, triggered, changed, triggered_step))
    any(view(changed, 1:N)) || return false
    forces_neighbours(mood.strategy) || return true
    copyto!(snapshot, 1, orders, 1, N)
    force_neighbours!(mood.strategy, pg, triggered, snapshot, orders, changed)
    return true
end

if COUNT
@eval HCC function (m::MOOD{<:MOODStrategy, MOODu2{T}})(
    g::DivergenceInterpolator, p_idx::Int, rho_i::State{M, T}, nb_slice::OrdinalRange{Int, Int},
    newRho::State{M, T}, pg::ParticleGrid{D, M, T}, U_prev::AbstractVector{State{M, T}}) where {D, M, T}
    δ = m.criterion.d
    Main.inc(:calls)
    nbs = neighbour_states(U_prev, pg)
    _dmp_fails(δ, rho_i, nb_slice, newRho, nbs) || return false
    Main.inc(:dmp_fail)
    if D == 2
        xp = pg.core.positions[p_idx]
        rbin = min(length(Main.RHIST), 1 + floor(Int, sqrt((xp[1] - 5)^2 + (xp[2] - 5)^2) / 0.25))   # first steps: centre ≈ (5, 5)
        Threads.atomic_add!(Main.RHIST[rbin], 1)
    end
    minU, maxU = findLocalExtrema(rho_i, nb_slice, nbs)
    viol = ntuple(mm -> abs(maxU[mm] - minU[mm]) >= δ^3 && (newRho[mm] < minU[mm] - δ || newRho[mm] > maxU[mm] + δ), Val(M))
    for mm in 1:M; viol[mm] && Threads.atomic_add!(Main.DMPCOMP[mm], 1); end
    if viol[1]                                     # density: overshoot and local range in units of δ (log10 bins from 1e-2)
        ex = max(newRho[1] - maxU[1], minU[1] - newRho[1])
        Threads.atomic_add!(Main.EXHIST[clamp(3 + floor(Int, log10(ex / δ)), 1, 12)], 1)
        Threads.atomic_add!(Main.RGHIST[clamp(3 + floor(Int, log10(max(maxU[1] - minU[1], 1e-300) / δ)), 1, 12)], 1)
    end
    if Main.VARIANT_PLATEAU && all(mm -> !viol[mm] || (maxU[mm] - minU[mm]) < Main.ETA * (Main.GHI[][mm] - Main.GLO[][mm]), 1:M)
        Main.inc(:accept); return false                  # every violating component is a plateau: no curvature needed
    end
    has_cached_gradients(g) || return true
    if _extract_order(g) < 3 || pg.core.particle_orders[p_idx] < 3
        Main.inc(:low_order); return true
    end
    grad_vec = cached_gradients(g); neighbors = pg.neighbor.indices
    first_fail = nothing; ok_on_violating = true
    h2 = minimum(pg.meta.dx)^2
    for d in 1:D
        c_i = grad_vec[p_idx][D + d]
        mini, maxi, minAbs, maxAbs = findLocalExtremaAbs(c_i, D + d, nb_slice, neighbors, grad_vec)
        for mm in 1:M
            Main.VARIANT_COMPS && !viol[mm] && continue          # (b) only the components that violated the DMP
            # (p) plateau: the component varies across the stencil by less than η of its global range
            Main.VARIANT_PLATEAU && (maxU[mm] - minU[mm]) < Main.ETA * (Main.GHI[][mm] - Main.GLO[][mm]) && continue
            ratio = maxAbs[mm] < T(1e-12) ? one(T) : minAbs[mm] / maxAbs[mm]
            if Main.VARIANT_UNITS                                 # (a) curvature as variation over one spacing, c h²
                sign_ok  = (mini[mm] * h2) * (maxi[mm] * h2) > -δ^2
                plateau  = maxAbs[mm] * h2 < δ
            else
                sign_ok  = mini[mm] * maxi[mm] > -δ
                plateau  = maxAbs[mm] < δ
            end
            ratio_ok = Main.VARIANT_NORATIO || ratio >= T(0.5) || plateau   # (c) no ratio test
            valid = sign_ok && ratio_ok
            Main.VARIANT_UNITS && plateau && (valid = true)
            if !valid
                viol[mm] && (ok_on_violating = false)
                first_fail === nothing && (first_fail = (viol[mm], sign_ok))
            end
        end
    end
    ok_on_violating && Main.inc(:would_accept_violating_only)
    if first_fail === nothing
        Main.inc(:accept); return false
    end
    v, s = first_fail
    Main.inc(s ? (v ? :ratio_violating : :ratio_other) : (v ? :sign_violating : :sign_other))
    return true
end
end

