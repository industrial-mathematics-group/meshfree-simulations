# Implosion (Hui, Li & Li 1999; Liska & Wendroff 2003): [0, 0.3]², slip walls on all sides, γ = 1.4, t = 2.5.
# (ρ, u, v, p) = (0.125, 0, 0, 0.14) for x + y < 0.15, (1, 0, 0, 1) elsewhere. The shock runs into the corner (0, 0),
# reflects, and a jet forms along the diagonal; the solution is symmetric under x ↔ y, which many schemes break.
# Lattice points (jitter 0), so the mirror image (y, x) of every particle is a particle again and the symmetry error
# is measured exactly by the post-processing statistic (implosion_symmetry!):
#   :symmetry_error        time series  Σ_i |ρ(x_i, y_i) − ρ(y_i, x_i)| / Σ_i ρ_i
#   :final_symmetry_error  scalar       the same at t_end
isdefined(@__MODULE__, :BENCH_METHODS) || Base.include(@__MODULE__, joinpath(@__DIR__, "bench_common.jl"))

"State `inner` where x + y < c, `outer` elsewhere."
struct DiagonalRiemannIC{T} <: HyperCloud.InitialCondition
    inner::State{4, T}
    outer::State{4, T}
    c::T
end
(ic::DiagonalRiemannIC{T})(pos) where {T} = pos[1] + pos[2] < ic.c ? ic.inner : ic.outer
function HyperCloud.build_ic(::Val{:diagonal_riemann}, conf::Dict, ctx::Dict)
    T = ctx[:T]
    return DiagonalRiemannIC{T}(State{4, T}(T.(conf[:inner])), State{4, T}(T.(conf[:outer])), T(conf[:c]))
end

function implosion_symmetry!(data::HyperCloud.AbstractSimData{D, DS, M, T}) where {D, DS, M, T}
    haskey(data.stats, :final_symmetry_error) && return false
    any(!iszero, data.params[:Grid_randomness_factor]) && return false     # needs lattice points
    h = Float64.(data.domain.spacing[1:DS])
    err = map(eachindex(data.t)) do n
        x, u = data.x[n], data.u[n]
        key(p) = (round(Int, p[1] / h[1]), round(Int, p[2] / h[2]))
        idx = Dict(key(x[i]) => i for i in eachindex(x))
        num = 0.0; den = 0.0
        for i in eachindex(x)
            j = get(idx, key((x[i][2], x[i][1])), 0)
            j == 0 && continue
            num += abs(Float64(u[i][1]) - Float64(u[j][1])); den += Float64(u[i][1])
        end
        T(num / den)
    end
    add_stat!(data, :symmetry_error, [SVector{M, T}(ntuple(_ -> e, M)) for e in err], :time)
    add_stat!(data, :final_symmetry_error, SVector{M, T}(ntuple(_ -> last(err), M)), Symbol[])
    return true
end

function get_config()
    γ = 1.4
    shared = bench_shared(:tmax => 2.5, :snapshots => 25, :Grid_mins => (0.0, 0.0), :Grid_maxs => (0.3, 0.3),
                          :Grid_Ns => (400, 400), :Grid_randomness_factor => (0.0, 0.0), :Grid_bc_map => bench_walls(),
                          :IC_name => :diagonal_riemann, :IC_inner => bench_cons(γ, 0.125, 0, 0, 0.14),
                          :IC_outer => bench_cons(γ, 1.0, 0, 0, 1.0), :IC_c => 0.15)
    return bench_config(shared; postprocess = :implosion_symmetry!, sources = [@__FILE__])
end
