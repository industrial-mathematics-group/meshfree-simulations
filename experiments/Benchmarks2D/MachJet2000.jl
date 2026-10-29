# MachJet2000: Mach 2000 jet (Zhang & Shu 2010): u_jet = 800, [0, 1] × [-0.25, 0.25], 640 × 320, t = 0.001.
isdefined(@__MODULE__, :JetInflowBC) || Base.include(@__MODULE__, joinpath(@__DIR__, "jet_common.jl"))
get_config() = jet_config(800.0, (0.0, -0.25), (1.0, 0.25), (640, 320), 0.001, @__FILE__)
