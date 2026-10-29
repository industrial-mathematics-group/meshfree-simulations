# GreshoLowMach: the Gresho vortex at peak Mach 0.01 (p0 = 1 / (γ 0.01²)), 128², t = 1 (0.8 rotations of the core; the
# time step follows the sound speed, 100× the flow speed), see gresho_common.jl; errors against the exact solution.
isdefined(@__MODULE__, :GreshoIC) || Base.include(@__MODULE__, joinpath(@__DIR__, "gresho_common.jl"))
get_config() = gresho_config(1 / (1.4 * 0.01^2), 1.0, (128, 128), @__FILE__)
