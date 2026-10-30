# MachJet80_Float64: the methods of MachJet80 whose GPU Float32 run ended with negative pressure, rerun in Float64 on the CPU
# (local, tools/local_queue.sh, queue cluster/queues/2026-10-30-bench2d-f64-local.txt); same parameters otherwise, so
# `@sim MachJet80_Float64` finds them. Report: reports/2026-10-30-benchmarks2d.md.
isdefined(@__MODULE__, :JetInflowBC) || Base.include(@__MODULE__, joinpath(@__DIR__, "jet_common.jl"))

function get_config()
    config = jet_config(30.0, (0.0, -0.5), (2.0, 0.5), (800, 400), 0.07, @__FILE__)
    config.shared[:T] = Float64
    config.defaults = Symbol.(["RK2MUSCL2(VK,AR,MLS2)"])
    return config
end
