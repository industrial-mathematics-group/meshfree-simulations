# Shared setup of the conservation-defect experiments (paper section "conservation defect": the floor of the shock
# test, frozen / coupled divergence degree, irregularity, support and weight, smooth periodic data).
# Included by the configs of this folder; each config is one paragraph of that section:
#
#   defect_floor     unperturbed lattice, every scheme                          (the floor of the shock test)
#   defect_shock     θ = 0.4, κ = 5.5: k_div = 2 / coupled (0) / K, EPD1/EPD2, aligned/geometric flux, limiters
#   defect_jitter    θ = 0, 0.1, …, 0.5 for the unstabilized K = 2, 4 and MUSCL2/4(EPD1)
#   defect_range     κ = 3.5 … 9.5 for MUSCL2(EPD1) and MUSCL4(EPD1, k_div = 4), aligned and geometric flux
#   defect_alpha     α = 0.5 … 4 for MUSCL4(EPD1, k_div = 4)
#   defect_smooth1d  periodic Gaussian advection (1-D), θ = 0.2, 0.4
#   defect_smooth2d  periodic 2-D Burgers before the shock (Table 7), θ = 0, 0.2, 0.4, κ = 3.5, 5.5
#
# Methods are those of experiments/methods_review.jl (SSP steppers, MOOD U2 with δ = 1e-4 Δx^D), labelled
# "<stepper>MUSCL<K>(<variant>,<AR|RU>,MLS<k_div>)": AR = aligned interface flux, RU = geometric (plain Rusanov),
# MLS 0 = divergence degree coupled to the MOOD order, MLS k = frozen at k. The defect is attached as a statistic
# by `conservation_defect!` (defect_stats.jl), see there for its definition.
isdefined(@__MODULE__, :add_review_methods!) ||
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "methods_review.jl"))

const DEFECT_SOURCES = [joinpath(@__DIR__, "defect_stats.jl")]
const DEFECT_N_1D = [ceil(Int, N) for N in 10 .^ (1.5:0.1:3.5)]        # 32 … 3163
const DEFECT_N_2D = [ceil(Int, N) for N in 10 .^ (1.5:0.1:2.3)]        # 32 … 200

# Method labels for (K, variant, flux, k_div), flux ∈ (:AR, :RU)
defect_label(K, variant, flux, mls) = review_label(K, variant, flux == :AR ? :AlignedRusanov : :Rusanov, mls)
defect_labels(specs) = [defect_label(s...) for s in specs]

function defect_methods()
    methods = get_master_method_dict()
    add_review_methods!(methods)
    return methods
end

# The 1-D Burgers shock: u_L = 1 → u_R = 0 at x = −2 on [−5, 5], outflow boundaries, t = 10 (the shock reaches
# x = 3), θ = 0.4, κ = 5.5, α = 1 (experiments/Conservation/conservation_burgers.jl).
defect_shock_shared() = Dict{Symbol, Any}(
    :PDE_name => :burgers,
    :PDE_D => 1,
    :T => Float64,
    :sim_func_name => :run_direct_simulation,
    :tmax => 10.0,
    :snapshots => 5,
    :Grid_domain => :rectangular,
    :Grid_mins => (-5.0,),
    :Grid_maxs => (5.0,),
    :Grid_Ns => (500,),
    :Grid_bc => Dict(0 => :outflow, 1 => :outflow),
    :Grid_periodic => false,
    :Grid_randomness_factor => (0.4,),
    :Grid_seed => 42,
    :Grid_mover => :none,
    :Grid_remove_ghosts => false,
    :Weight_name => :exponential,
    :Weight_alpha => 1.0,
    :Weight_range => 5.5,
    :IC_name => :riemann,
    :IC_uL => (1.0,),
    :IC_uR => (0.0,),
    :IC_p0 => (-2.0,),
    :IC_n0 => (1.0,),
    :Time_CFL => 0.2,
)

# Builds the config of one paragraph: `specs` are (K, variant, flux, k_div) tuples, `varied` extra sweeps
# (the grid sweep is added).
function defect_config(shared, specs; varied = Dict{Symbol, Any}(), Ns = DEFECT_N_1D)
    v = create_varied_dict()
    push!(v, :Grid_Ns__all => Ns)
    for (k, vals) in varied
        push!(v, k => vals)
    end
    return SimulationConfig(shared, defect_methods(), defect_labels(specs); varied = v,
                            reference = "analytical_solution", postprocess = :conservation_defect!,
                            sources = DEFECT_SOURCES)
end
