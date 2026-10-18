# Standalone reimplementation of the library's MUSCL/Rusanov scheme on a HyperCloudCore point cloud
# (scalar laws), with alternative divergence formulations:
#
#   :gfd   library form   dU_i/dt = -2 Σ_j c_ij · (F*_ij - F_i)           (consistent, not conservative)
#   :sym   face form      dU_i/dt = -(1/V_i) Σ_j A_ij · F*_ij,
#                         A_ij = V_i c_ij - V_j c_ji                       (conservative, Σ_j A_ij ≠ 0)
#   :symc  :sym with A_ij + (μ_j - μ_i) such that Σ_j A_ij = 0 for every i
#          (one graph-Laplacian solve per direction and geometry)          (conservative, free-stream exact)
#
# c_ij are the MLS first-derivative coefficients of the library fit (degree `deg`, no constant term).
# Σ_i V_i U_i is conserved exactly by :sym and :symc for the volumes V used to build A.

using SparseArrays

struct Pairs{D, T}
    n::Int
    ranges::Vector{UnitRange{Int}}
    J::Vector{Int}                 # neighbour of slot k
    dist::Vector{SVector{D, T}}    # x_j - x_i (periodic)
    C::Vector{SVector{D, T}}       # divergence coefficients c_ij (degree `deg`)
    G::Vector{SVector{D, T}}       # gradient coefficients (degree 1, MUSCL2 reconstruction)
    rev::Vector{Int}               # slot of (j, i)
    V::Vector{T}
end

function Pairs(pg::ParticleGrid{D, M, T}; deg = 1, V = nothing) where {D, M, T}
    nb = pg.neighbor; n = pg.meta.N
    ranges = [nb.ranges[i] for i in 1:n]
    J = collect(nb.indices[1:last(ranges[n])]); dist = collect(nb.distances[1:last(ranges[n])])
    _, _, C = coefficients(pg, deg)
    _, _, G = coefficients(pg, 1)
    slot = Dict{Tuple{Int, Int}, Int}()
    for i in 1:n, k in ranges[i]; slot[(i, J[k])] = k; end
    rev = [slot[(J[k], i)] for i in 1:n for k in ranges[i]]
    VV = V === nothing ? collect(pg.core.volumes[1:n]) : V
    return Pairs{D, T}(n, ranges, J, dist, C, G, rev, VV)
end

"Face coefficients A_ij for `form` (:sym, :symc, :symc1); nothing for :gfd."
function face_coefficients(P::Pairs{D, T}, form::Symbol) where {D, T}
    form === :gfd && return nothing
    form === :symc1 && return face_coefficients_exact(P)
    form === :symc1u && return face_coefficients_exact(P; weighted = false)
    A = [P.V[i] * P.C[k] - P.V[P.J[k]] * P.C[P.rev[k]] for i in 1:P.n for k in P.ranges[i]]
    form === :sym && return A
    # Minimal antisymmetric correction δ_ij = μ_j - μ_i with Σ_j (A_ij + δ_ij) = 0:  (L μ)_i = r_i
    rows = Int[]; cols = Int[]; vals = T[]
    for i in 1:P.n
        push!(rows, i); push!(cols, i); push!(vals, T(length(P.ranges[i])))
        for k in P.ranges[i]; push!(rows, i); push!(cols, P.J[k]); push!(vals, -one(T)); end
    end
    L = sparse(rows, cols, vals, P.n, P.n)
    L[1, :] .= 0; L[1, 1] = 1                      # fix the constant (L is singular, Σ r = 0)
    F = lu(L)
    r = [sum(A[k] for k in P.ranges[i]) for i in 1:P.n]
    μ = [F \ [i == 1 ? zero(T) : r[i][d] for i in 1:P.n] for d in 1:D]
    return [A[k] + SVector{D, T}(ntuple(d -> μ[d][P.J[k]] - μ[d][i], Val(D))) for i in 1:P.n for k in P.ranges[i]]
end

"""
    face_coefficients_exact(P) -> A

Antisymmetric face coefficients with zero row sums AND linear exactness
`Σ_j A_ij ⊗ (x_j - x_i) = 2 V_i I`: the minimal-norm antisymmetric correction of
`V_i c_ij - V_j c_ji` satisfying these (D + D²) constraints per particle (one sparse solve).
"""
function face_coefficients_exact(P::Pairs{D, T}; weighted = true) where {D, T}
    A0 = [P.V[i] * P.C[k] - P.V[P.J[k]] * P.C[P.rev[k]] for i in 1:P.n for k in P.ranges[i]]
    # unordered pairs p = (i < j): unknown δ_p ∈ R^D, δ_ij = δ_p, δ_ji = -δ_p
    pid = zeros(Int, length(P.J)); np = 0
    for i in 1:P.n, k in P.ranges[i]
        j = P.J[k]
        if i < j; np += 1; pid[k] = np; pid[P.rev[k]] = -np; end
    end
    nc = D + D * D                                       # constraints per particle
    rows = Int[]; cols = Int[]; vals = T[]; rhs = zeros(T, nc * P.n)
    for i in 1:P.n
        r = zero(SVector{D, T}); Mo = zeros(T, D, D)
        for k in P.ranges[i]; r += A0[k]; Mo .+= A0[k] * P.dist[k]'; end
        E = 2 * P.V[i] * Matrix{T}(LinearAlgebra.I, D, D) .- Mo
        base = nc * (i - 1)
        for d in 1:D; rhs[base + d] = -r[d]; end
        for d in 1:D, e in 1:D; rhs[base + D + (d - 1) * D + e] = E[d, e]; end
        for k in P.ranges[i]
            p = abs(pid[k]); sgn = sign(pid[k]); x = P.dist[k]
            for d in 1:D
                col = D * (p - 1) + d
                push!(rows, base + d); push!(cols, col); push!(vals, sgn)
                for e in 1:D
                    push!(rows, base + D + (d - 1) * D + e); push!(cols, col); push!(vals, sgn * x[e])
                end
            end
        end
    end
    B = sparse(rows, cols, vals, nc * P.n, D * np)
    # weighted minimal norm Σ δ_p² / ω_p, ω_p = |A0_p|²: small coefficients get small corrections
    ω = ones(T, D * np)
    if weighted
        for i in 1:P.n, k in P.ranges[i]
            pid[k] > 0 || continue
            for d in 1:D; ω[D * (pid[k] - 1) + d] = sum(abs2, A0[k]) + eps(T); end
        end
    end
    W = spdiagm(ω)
    BWBt = B * W * B'
    λ = (BWBt + T(1e-12) * maximum(abs, BWBt) * sparse(LinearAlgebra.I, size(BWBt)...)) \ rhs
    δ = W * (B' * λ)
    res = norm(B * δ - rhs) / norm(rhs)
    res > 1e-6 && @warn "exact face coefficients: constraint residual $res"
    return [A0[k] + sign(pid[k]) * SVector{D, T}(ntuple(d -> δ[D * (abs(pid[k]) - 1) + d], Val(D))) for i in 1:P.n for k in P.ranges[i]]
end

# Rusanov flux along the unit normal n (states: i side, j side), scalar laws: A·F* = |A| f̂(n)
@inline function rusanov_normal(eq, ui::T, uj::T, A::SVector{D, T}) where {D, T}
    a = norm(A); a == 0 && return zero(T)
    n = A / a
    Fi = flux(eq, State{1, T}(ui)); Fj = flux(eq, State{1, T}(uj))
    Fn_i = sum(n[d] * Fi[d][1] for d in 1:D); Fn_j = sum(n[d] * Fj[d][1] for d in 1:D)
    s = max(sum(abs(n[d]) * max_eigenvalue(eq, State{1, T}(ui), d) for d in 1:D),
            sum(abs(n[d]) * max_eigenvalue(eq, State{1, T}(uj), d) for d in 1:D))
    return a * ((Fn_i + Fn_j) - s * (uj - ui)) / 2
end

"Scale of the Rusanov dissipation (1 = library; 0 = central flux, for the order-reduction test)."
const DISSIPATION = Ref(1.0)

# Rusanov interface flux (per direction, states sorted by the sign of the distance as in the library)
@inline function rusanov(eq, ui::T, uj::T, d::SVector{D, T}) where {D, T}
    Fi = flux(eq, State{1, T}(ui)); Fj = flux(eq, State{1, T}(uj))
    SVector{D, T}(ntuple(Val(D)) do a
        uL, uR = d[a] > 0 ? (ui, uj) : (uj, ui)
        FL, FR = d[a] > 0 ? (Fi[a][1], Fj[a][1]) : (Fj[a][1], Fi[a][1])
        s = max(max_eigenvalue(eq, State{1, T}(uL), a), max_eigenvalue(eq, State{1, T}(uR), a))
        (FL + FR - DISSIPATION[] * s * (uR - uL)) / 2
    end)
end

"""
    rhs!(K, P, A, eq, u; order) -> K

Semi-discrete right-hand side. `order = 1` (piecewise constant) or `2` (unlimited MUSCL with
degree-1 MLS gradients), `A` from `face_coefficients` (`nothing` = library form).
"""
function rhs!(K, P::Pairs{D, T}, A, eq, u; order = 2, limiter = :none, riemann = :axis) where {D, T}
    n = P.n
    g = Vector{SVector{D, T}}(undef, n)
    for i in 1:n
        acc = zero(SVector{D, T})
        if order == 2
            for k in P.ranges[i]; acc += P.G[k] * (u[P.J[k]] - u[i]); end
            if limiter === :bj                      # Barth-Jespersen on the interface values
                umax = u[i]; umin = u[i]
                for k in P.ranges[i]; umax = max(umax, u[P.J[k]]); umin = min(umin, u[P.J[k]]); end
                phi = one(T)
                for k in P.ranges[i]
                    du = dot(acc, P.dist[k]) / 2
                    du > 0 && (phi = min(phi, (umax - u[i]) / du))
                    du < 0 && (phi = min(phi, (umin - u[i]) / du))
                end
                acc *= max(phi, zero(T))
            end
        end
        g[i] = acc
    end
    for i in 1:n
        Fi = SVector{D, T}(ntuple(a -> flux(eq, State{1, T}(u[i]))[a][1], Val(D)))
        acc = zero(T)
        for k in P.ranges[i]
            j = P.J[k]; d = P.dist[k]
            ui = u[i] + dot(g[i], d) / 2; uj = u[j] - dot(g[j], d) / 2
            if A !== nothing && riemann === :normal
                acc += rusanov_normal(eq, ui, uj, A[k]) / P.V[i]
            else
                Fs = rusanov(eq, ui, uj, d)
                acc += A === nothing ? 2 * dot(P.C[k], Fs - Fi) : dot(A[k], Fs) / P.V[i]
            end
        end
        K[i] = -acc
    end
    return K
end

"SSP-RK2 to `tmax` with the library CFL estimate; returns (u, steps)."
function solve(P::Pairs{D, T}, A, eq, pg, u0, tmax; cfl = 0.2, order = 2, limiter = :none, riemann = :axis, interp = HyperCloud.NoFallbackGrad()) where {D, T}
    u = copy(u0); K1 = similar(u); K2 = similar(u); t = zero(T); k = 0
    while t < tmax
        copyto!(pg.rhos, 1, State{1, T}.(u), 1, P.n)
        dt = min(cfl * get_time_step(pg, eq, interp), tmax - t)
        rhs!(K1, P, A, eq, u; order, limiter, riemann); u1 = u .+ dt .* K1
        rhs!(K2, P, A, eq, u1; order, limiter, riemann); u .= (u .+ u1 .+ dt .* K2) ./ 2
        t += dt; k += 1
    end
    return u, k
end

# -------------------------------------------------------------------------------------------------
# Exact smooth solutions and a fine 1-D reference
# -------------------------------------------------------------------------------------------------
"Burgers with flux u²/2 in every direction: u = u0(x - u t (1,…,1)), solved by Newton (pre-shock)."
function burgers_exact(u0, x::SVector{D, T}, t) where {D, T}
    u = u0(x)
    for _ in 1:50
        y = x .- u * t
        h = 1e-7
        f = u - u0(y); fp = 1 + t * sum((u0(y .+ h * SVector{D, T}(ntuple(b -> b == a, Val(D)))) - u0(y)) / h for a in 1:D)
        du = f / fp; u -= du
        abs(du) < 1e-14 && break
    end
    return u
end

"Godunov finite volumes for 1-D Burgers on [0,1) (periodic, N cells, CFL 0.4): reference solution."
function burgers_fv_1d(u0, N, tmax)
    dx = 1 / N; x = ((1:N) .- 0.5) .* dx
    u = [sum(u0(SVector(xi - dx / 2 + dx * (q + 0.5) / 8)) for q in 0:7) / 8 for xi in x]
    god(l, r) = l <= r ? (l > 0 ? l^2 / 2 : r < 0 ? r^2 / 2 : 0.0) : max(l^2, r^2) / 2
    t = 0.0
    while t < tmax
        dt = min(0.4 * dx / maximum(abs, u), tmax - t)
        F = [god(u[i], u[mod1(i + 1, N)]) for i in 1:N]
        u = u .- dt / dx .* (F .- F[[N; 1:N-1]])
        t += dt
    end
    return x, u
end

# =================================================================================================
# Higher order: degree-q MLS reconstruction, degree-q GFD divergence, and conservative face
# coefficients exact for flux polynomials of degree p at the midpoints
# =================================================================================================
struct PairsHO{D, T}
    n::Int
    ranges::Vector{UnitRange{Int}}
    J::Vector{Int}
    dist::Vector{SVector{D, T}}
    S::Vector{Vector{T}}           # scaled MLS coefficients of slot k (degree q, basis without constant)
    rev::Vector{Int}
    V::Vector{T}
    q::Int
    invL::T
end

function PairsHO(pg::ParticleGrid{D, M, T}; q = 1) where {D, M, T}
    nb = pg.neighbor; n = pg.meta.N; invL = one(T) / minimum(pg.meta.dx)
    ranges = [nb.ranges[i] for i in 1:n]
    J = collect(nb.indices[1:last(ranges[n])]); dist = collect(nb.distances[1:last(ranges[n])])
    S = Vector{Vector{T}}(undef, length(J))
    for i in 1:n
        Minv, _ = _inv_moment(Val(q), ranges[i], nb.distances, nb.weights, invL, Val(D), T)
        for k in ranges[i]
            S[k] = Minv * collect(HCC.build_basis(Val(q), dist[k] * invL)) * nb.weights[k]
        end
    end
    slot = Dict{Tuple{Int, Int}, Int}()
    for i in 1:n, k in ranges[i]; slot[(i, J[k])] = k; end
    rev = [slot[(J[k], i)] for i in 1:n for k in ranges[i]]
    return PairsHO{D, T}(n, ranges, J, dist, S, rev, collect(pg.core.volumes[1:n]), q, invL)
end

"GFD divergence coefficients c_ij (physical) of the degree-q fit."
gfd_coefficients(P::PairsHO{D, T}) where {D, T} = [SVector{D, T}(ntuple(d -> P.S[k][d] * P.invL, Val(D))) for k in eachindex(P.J)]

"""
    face_coefficients_p(P::PairsHO, p) -> A

Minimal-norm antisymmetric correction of `V_i c_ij - V_j c_ji` such that
`(1/V_i) Σ_j A_ij · F(x_i + d_ij/2)` is exact for flux polynomials of degree ≤ p.
"""
function face_coefficients_p(P::PairsHO{D, T}, p::Int) where {D, T}
    C = gfd_coefficients(P)
    A0 = [P.V[i] * C[k] - P.V[P.J[k]] * C[P.rev[k]] for i in 1:P.n for k in P.ranges[i]]
    pid = zeros(Int, length(P.J)); np = 0
    for i in 1:P.n, k in P.ranges[i]
        j = P.J[k]
        if i < j; np += 1; pid[k] = np; pid[P.rev[k]] = -np; end
    end
    B = length(HCC.build_basis(Val(p), zero(SVector{D, T})))
    nc = D * (1 + B)                                     # per particle: (moment 0 and moments b) × direction e
    rows = Int[]; cols = Int[]; vals = T[]; rhs = zeros(T, nc * P.n)
    for i in 1:P.n
        base = nc * (i - 1)
        cur = zeros(T, nc)
        for k in P.ranges[i]
            m = vcat(one(T), collect(HCC.build_basis(Val(p), P.dist[k] / 2)))
            for e in 1:D, a in 1:B + 1; cur[(e - 1) * (B + 1) + a] += A0[k][e] * m[a]; end
            pp = abs(pid[k]); sgn = sign(pid[k])
            for e in 1:D, a in 1:B + 1
                push!(rows, base + (e - 1) * (B + 1) + a); push!(cols, D * (pp - 1) + e); push!(vals, sgn * m[a])
            end
        end
        for e in 1:D, a in 1:B + 1
            target = (a - 1 == e) ? P.V[i] : zero(T)     # moment b = e (the linear term along e) equals V_i
            rhs[base + (e - 1) * (B + 1) + a] = target - cur[(e - 1) * (B + 1) + a]
        end
    end
    Bm = sparse(rows, cols, vals, nc * P.n, D * np)
    ω = ones(T, D * np)
    for i in 1:P.n, k in P.ranges[i]
        pid[k] > 0 || continue
        for d in 1:D; ω[D * (pid[k] - 1) + d] = sum(abs2, A0[k]) + eps(T); end
    end
    W = spdiagm(ω); BWBt = Bm * W * Bm'
    λ = (BWBt + T(1e-12) * maximum(abs, BWBt) * sparse(LinearAlgebra.I, size(BWBt)...)) \ rhs
    δ = W * (Bm' * λ)
    res = norm(Bm * δ - rhs) / norm(rhs)
    res > 1e-6 && @warn "face coefficients (p = $p): constraint residual $res"
    return [A0[k] + sign(pid[k]) * SVector{D, T}(ntuple(d -> δ[D * (abs(pid[k]) - 1) + d], Val(D))) for i in 1:P.n for k in P.ranges[i]]
end

"RHS with degree-q reconstruction to the midpoints; A = nothing: GFD form with c of the same fit."
function rhs_ho!(K, P::PairsHO{D, T}, A, C, eq, u; riemann = :axis) where {D, T}
    n = P.n; Bq = length(P.S[1])
    a = [zeros(T, Bq) for _ in 1:n]
    for i in 1:n, k in P.ranges[i]
        a[i] .+= P.S[k] .* (u[P.J[k]] - u[i])
    end
    for i in 1:n
        Fi = SVector{D, T}(ntuple(e -> flux(eq, State{1, T}(u[i]))[e][1], Val(D)))
        acc = zero(T)
        for k in P.ranges[i]
            j = P.J[k]; d = P.dist[k]
            bi = HCC.build_basis(Val(P.q), d * (P.invL / 2)); bj = HCC.build_basis(Val(P.q), -d * (P.invL / 2))
            ui = u[i] + dot(a[i], bi); uj = u[j] + dot(a[j], bj)
            if A === nothing
                acc += 2 * dot(C[k], rusanov(eq, ui, uj, d) - Fi)
            elseif riemann === :normal
                acc += rusanov_normal(eq, ui, uj, A[k]) / P.V[i]
            else
                acc += dot(A[k], rusanov(eq, ui, uj, d)) / P.V[i]
            end
        end
        K[i] = -acc
    end
    return K
end

"Classical RK4 (smooth tests) with the library CFL estimate of degree `cfl_deg`."
function solve_ho(P::PairsHO{D, T}, A, eq, pg, u0, tmax; cfl = 0.2, riemann = :normal, stepper = :rk4,
                  C = gfd_coefficients(P)) where {D, T}
    u = copy(u0); t = zero(T); k = 0
    K = [similar(u) for _ in 1:4]
    interp = MUSCL(T, D, 1, P.q + 1, NoLimiter(), RusanovFlux())
    f!(Kx, ux) = rhs_ho!(Kx, P, A, C, eq, ux; riemann)
    while t < tmax
        copyto!(pg.rhos, 1, State{1, T}.(u), 1, P.n)
        dt = min(cfl * get_time_step(pg, eq, interp), tmax - t)
        if stepper === :rk4
            f!(K[1], u); f!(K[2], u .+ dt / 2 .* K[1]); f!(K[3], u .+ dt / 2 .* K[2]); f!(K[4], u .+ dt .* K[3])
            u .+= dt / 6 .* (K[1] .+ 2 .* K[2] .+ 2 .* K[3] .+ K[4])
        else
            f!(K[1], u); u1 = u .+ dt .* K[1]; f!(K[2], u1); u .= (u .+ u1 .+ dt .* K[2]) ./ 2
        end
        t += dt; k += 1
    end
    return u, k
end
