isdefined(@__MODULE__, :defect_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "defect_setup.jl"))
# The floor of the shock test: on the unperturbed lattice every scheme gives the same defect, zero to round-off or
# −h u_L / (2 m(t_end)) depending on N (half a particle of the inflow state).
function get_config()
    shared = defect_shock_shared()
    shared[:Grid_randomness_factor] = (0.0,)
    specs = [[(K, :none, :AR, K) for K in 2:5];
             [(K, v, :AR, 2) for K in 2:5 for v in (:EPD1, :EPD2)];
             [(K, :EPD1, :AR, K) for K in 2:5]]
    return defect_config(shared, unique(specs))
end
