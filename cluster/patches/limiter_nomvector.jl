# Test patch for HyperCloudCore 0.6.2 (not part of the package): Barth-Jespersen / Venkatakrishnan slope
# limiting without MVectors. In the original, the final `ntuple` closure captures the MVector `phi_i`, so it
# is heap allocated; on CUDA that is a device malloc per call, which exhausts the device heap
# ("Out of dynamic GPU memory (trying to allocate 16 bytes)"). Only immutable SVectors here, and the
# closures only capture variables that are assigned once (others would be boxed).
@eval HyperCloud.HyperCloudCore function _limit_slopes(
    strategy::Union{BarthJespersenLimiter{Mode}, VenkatakrishnanLimiter{Mode}},
    raw_grad::SVector{B_LEN, State{M, T}},
    nb_slice::OrdinalRange{Int, Int},
    f_i::State{M, T},
    f_neighbors::AbstractVector{State{M, T}},
    pg::ParticleGrid{D, M, T},
    distVec::AbstractVector{Space{D, T}},
    ::Val{DEGREE}
) where {B_LEN, M, D, T, DEGREE, Mode}
    u_max = f_i
    u_min = f_i
    @inbounds for global_idx in nb_slice
        f_j = f_neighbors[global_idx]
        u_max = max.(u_max, f_j)
        u_min = min.(u_min, f_j)
    end
    dmax = u_max - f_i
    dmin = u_min - f_i

    phi = ones(State{M, T})
    @inbounds for global_idx in nb_slice
        p_interface = build_basis(Val(DEGREE), T(0.5) * distVec[global_idx])
        delta_recon = zero(State{M, T})
        for k in 1:B_LEN
            delta_recon += raw_grad[k] * p_interface[k]
        end
        phi = map(phi, delta_recon, dmax, dmin) do ph, rec, dmx, dmn
            abs(rec) > T(1e-12) ? min(ph, limiter_phi(strategy, (rec > zero(T) ? dmx : dmn) / rec)) : ph
        end
    end

    phi_final = phi
    return SVector{B_LEN, State{M, T}}(ntuple(Val(B_LEN)) do k
        (Mode === :hard || k <= D) ? raw_grad[k] .* phi_final : raw_grad[k]
    end)
end
