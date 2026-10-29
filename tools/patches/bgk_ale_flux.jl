# Test patch (not part of HyperCloud): on a moving grid the scheme sees ALEFrame(eq, w), which forwards flux and wave
# speeds but not BGK's own RusanovFlux method (per-velocity upwind dissipation |v_k,d|). The generic Rusanov then
# dissipates every velocity with the largest relative speed (v_max + |w|). This method restores the per-velocity
# upwind dissipation in the moving frame, |v_k,d − w_d|. Loaded by tools/adapt_study.jl with ADAPT_BGK_ALE_FLUX=1.
@inline function (::HyperCloud.HyperCloudCore.RusanovFlux)(f_L::HyperCloud.HyperCloudCore.Flux{D, M, T}, f_R::HyperCloud.HyperCloudCore.Flux{D, M, T},
        F_L::HyperCloud.HyperCloudCore.Flux{D, M, T}, F_R::HyperCloud.HyperCloudCore.Flux{D, M, T},
        a::HyperCloud.HyperCloudCore.ALEFrame{D, M, T, <:HyperCloud.BGKEquation}) where {D, M, T}
    return HyperCloud._per_axis(a.eq) do v, d
        T(0.5) * (F_L[d] + F_R[d] - abs.(v .- a.w[d]) .* (f_R[d] - f_L[d]))
    end
end
