# MachJet80: Mach 80 jet (Ha et al. 2005, Zhang & Shu 2010): u_jet = 30, [0, 2] × [-0.5, 0.5], 800 × 400, t = 0.07.
isdefined(@__MODULE__, :JetInflowBC) || Base.include(@__MODULE__, joinpath(@__DIR__, "jet_common.jl"))
get_config() = jet_config(30.0, (0.0, -0.5), (2.0, 0.5), (800, 400), 0.07, @__FILE__)
