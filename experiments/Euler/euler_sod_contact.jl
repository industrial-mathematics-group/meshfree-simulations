# euler_sod_contact: the Euler counterpart of experiments/BGK/bgk_sod_contact.jl (study/bgk-lagrangian-contact). HyperCloud's
# euler_1d_sod (the same Sod setup as bgk_1d_sod: [-0.5, 0.5], jitter 0.2, seed 10, R = 5.5, CFL 0.2, t = 0.2), Float64,
# with the methods of bgk_sod.jl on static and Lagrangian particles (MLS 2 volumes), to see whether the contact smearing of
# ARS222MUSCL2(EPD2) on Lagrangian particles needs the kinetic model. Error statistics against the exact Riemann solution.
# Run locally: SIM_BACKEND=cpu SIM_T=Float64 SIM_SAVE=data/contact julia --project=env_web -t 12 tools/run_sim.jl euler_sod_contact
module _euler_sod_contact_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "BGK", "bgk_sod.jl"))
end

function get_config()
    shared = copy(euler_1d_sod().shared)
    shared[:T] = Float64
    shared[:snapshots] = 20
    base = _euler_sod_contact_Base.get_config().methods
    lag = Dict{Symbol, Any}(:Grid_mover => :lagrangian, :Grid_volumes => :mls, :Grid_volumes_degree => 2,
                            :Adapt_enabled => false)
    epd2 = merge(base[Symbol("ARS222MUSCL2(EPD2)")], lag)
    vk   = merge(base[Symbol("ARS222MUSCL2(VK)")], lag)
    variants = [
        "EPD2"              => epd2,
        "EPD2 static"       => merge(epd2, Dict{Symbol, Any}(:Grid_mover => :none)),
        "VK"                => vk,
        "VK static"         => merge(vk, Dict{Symbol, Any}(:Grid_mover => :none)),
        "EPD2 Rusanov"      => merge(epd2, Dict{Symbol, Any}(:Flux_name => :Rusanov)),
        "VK AlignedRusanov" => merge(vk, Dict{Symbol, Any}(:Flux_name => :AlignedRusanov)),
        "EPD1"              => merge(epd2, Dict{Symbol, Any}(:MOOD_strategy => :EPD1)),
        "EPD2 ARS233"       => merge(epd2, Dict{Symbol, Any}(:Time_stepper => :ARS233)),
        "EPD2 K3"           => merge(epd2, Dict{Symbol, Any}(:Scheme_order => 3)),
    ]
    methods = Dict{Symbol, Dict{Symbol, Any}}(Symbol(n) => m for (n, m) in variants)
    varied = create_varied_dict()
    push!(varied, :Grid_Ns__1 => [200, 400, 800, 1600, 3200])
    return SimulationConfig(shared, methods, [Symbol(n) for (n, _) in variants]; varied = varied,
                            reference = :analytical_solution, sources = [@__FILE__])
end
