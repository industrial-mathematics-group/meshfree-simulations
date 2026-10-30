# Noh_Float64: the methods of Noh whose GPU Float32 run ended with negative pressure, rerun in Float64 on the CPU
# (local, tools/local_queue.sh, queue cluster/queues/2026-10-30-bench2d-f64-local.txt); same parameters otherwise, so
# `@sim Noh_Float64` finds them. Report: reports/2026-10-30-benchmarks2d.md.
isdefined(@__MODULE__, :noh_config) || Base.include(@__MODULE__, joinpath(@__DIR__, "Noh.jl"))

function get_config()
    config = noh_config(@__FILE__)
    config.shared[:T] = Float64
    config.defaults = Symbol.(["RK2MUSCL2(VK,AR,MLS2)", "RK3MUSCL3(EPD1,AR,MLS3)", "RK3MUSCL3(EPD2,AR,MLS3)", "RK4SSPMUSCL4(EPD2,AR,MLS2)", "RK4SSPMUSCL4(EPD1,AR,MLS4)", "RK4SSPMUSCL4(EPD2,AR,MLS4)"])
    return config
end
