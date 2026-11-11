# bgk_sod_contact: one-factor tests of the contact smearing of ARS222MUSCL2(EPD2) on Lagrangian particles
# (study/bgk-lagrangian-contact). BGK Sod as in bgk_sod_sweep (τ = 1e-5, Float64, MLS 2 volumes, no particle
# management); every method is the Lagrangian EPD2 run ("EPD2") or the Lagrangian VK run ("VK") with one ingredient
# changed, named after the change.
# Run locally: SIM_BACKEND=cpu SIM_T=Float64 SIM_SAVE=data/contact julia --project=env_web -t 12 tools/run_sim.jl bgk_sod_contact
module _bgk_sod_contact_Base
    using HyperCloud, StaticArrays
    Base.include(@__MODULE__, joinpath(@__DIR__, "bgk_sod.jl"))
end

function get_config()
    config = _bgk_sod_contact_Base.get_config()
    config.shared[:T] = Float64
    config.shared[:PDE_tau] = 1e-5
    empty!(config.varied)
    config.varied[:Grid_Ns__1] = [800]
    lag = Dict{Symbol, Any}(:Grid_mover => :lagrangian, :Grid_volumes => :mls, :Grid_volumes_degree => 2,
                            :Adapt_enabled => false)
    epd2 = merge(config.methods[Symbol("ARS222MUSCL2(EPD2)")], lag)
    vk   = merge(config.methods[Symbol("ARS222MUSCL2(VK)")], lag)
    variants = [
        # references
        "EPD2"              => epd2,
        "EPD2 static"       => merge(epd2, Dict{Symbol, Any}(:Grid_mover => :none)),
        "VK"                => vk,
        "VK static"         => merge(vk, Dict{Symbol, Any}(:Grid_mover => :none)),
        # numerical flux (the VK method uses Rusanov, the EPD2 methods AlignedRusanov)
        "EPD2 Rusanov"      => merge(epd2, Dict{Symbol, Any}(:Flux_name => :Rusanov)),
        "VK AlignedRusanov" => merge(vk, Dict{Symbol, Any}(:Flux_name => :AlignedRusanov)),
        # MOOD: strategy and criterion
        "EPD1"              => merge(epd2, Dict{Symbol, Any}(:MOOD_strategy => :EPD1)),
        "EPD0"              => merge(epd2, Dict{Symbol, Any}(:MOOD_strategy => :EPD0)),
        "SEPD0"             => merge(epd2, Dict{Symbol, Any}(:MOOD_strategy => :SEPD0)),
        "EPD2 U1"           => merge(epd2, Dict{Symbol, Any}(:MOOD_criterion => :U1)),
        "EPD2 delta 1e-2"   => merge(epd2, Dict{Symbol, Any}(:MOOD_delta_relax => 1e-2)),
        "EPD2 plateau 1e-2" => merge(epd2, Dict{Symbol, Any}(:MOOD_plateau => 1e-2)),
        "no MOOD"           => merge(epd2, Dict{Symbol, Any}(:MOOD_criterion => :none)),
        # stencil radius (absolute, R times the initial spacing)
        "EPD2 R 2.5"        => merge(epd2, Dict{Symbol, Any}(:Weight_range => 2.5)),
        "EPD2 R 3.5"        => merge(epd2, Dict{Symbol, Any}(:Weight_range => 3.5)),
        "EPD2 R 7.5"        => merge(epd2, Dict{Symbol, Any}(:Weight_range => 7.5)),
        "VK R 3.5"          => merge(vk, Dict{Symbol, Any}(:Weight_range => 3.5)),
        # time stepper vs scheme order (ARS233MUSCL3(EPD2) was much better on the same grid)
        "EPD2 ARS233"       => merge(epd2, Dict{Symbol, Any}(:Time_stepper => :ARS233)),
        "EPD2 K3"           => merge(epd2, Dict{Symbol, Any}(:Scheme_order => 3)),
        "EPD2 CFL 0.1"      => merge(epd2, Dict{Symbol, Any}(:Time_CFL => 0.1)),
        # particles: no initial jitter, half-Lagrangian motion, Shepard volumes
        "EPD2 jitter 0"     => merge(epd2, Dict{Symbol, Any}(:Grid_randomness_factor => (0.0,))),
        "VK jitter 0"       => merge(vk, Dict{Symbol, Any}(:Grid_randomness_factor => (0.0,))),
        "EPD2 alpha 0.5"    => merge(epd2, Dict{Symbol, Any}(:Grid_mover => :scaled_lagrangian, :Grid_mover_alpha => 0.5)),
        "EPD2 shepard"      => merge(epd2, Dict{Symbol, Any}(:Grid_volumes => :shepard)),
    ]
    config.defaults = Symbol[]
    for (name, m) in variants
        config.methods[Symbol(name)] = m
        push!(config.defaults, Symbol(name))
    end
    push!(config.sources, @__FILE__)
    return config
end
