# E4: why ARS222MUSCL2(EPD2) on the Lagrangian grid loses the plateau between contact and shock (N = 800, τ = 1e-5).
#  a) when: the solution over time, Lagrangian vs static
#  b) where: the conserved components ρ, ρu, E at t = 0.2
#  c) the particles: spacing of the Lagrangian particles, EPD2 vs VK (same mover), from studio_data
include(joinpath(@__DIR__, "session.jl"))
using PDEStudio.CairoMakie   # PDEStudio's CairoMakie (not a direct dependency of env_web)
pair = [Symbol("ARS222MUSCL2(EPD2) static"), Symbol("ARS222MUSCL2(EPD2) lagrangian")]
studio_set!(:methods => pair)
studio_set!(:Grid_Ns__1 => 800, :PDE_tau => 1e-5, :compare_target => :time, :compare_columns => 4, :plot_style => :particles_lines,
            :plot_width => 360, :plot_height => 230)
studio_option!(:labels, :x_label, "x"); studio_option!(:labels, :y_label, "ρ")
studio_press!(:update_plot)
println(shot("e4a_time_epd2_static_vs_lagrangian"))

studio_set!(:compare_target => :component, :compare_link => :decoupled, :compare_columns => 3, :t => 0.2, :plot_width => 470, :plot_height => 280)
studio_option!(:labels, :y_label, "")
studio_press!(:update_plot)
println(shot("e4b_components_t0.2"))

# c) particle spacing (the Lagrangian particles move with the flow; MUSCL2 EPD2 vs VK)
studio_set!(:methods => [Symbol("ARS222MUSCL2(EPD2) lagrangian"), Symbol("ARS222MUSCL2(VK) lagrangian")])
studio_set!(:compare_target => :methods, :compare_columns => 2, :component => 1, :t => 0.2)
fig = Figure(size = (1100, 380), fontsize = 15)
ax1 = Axis(fig[1, 1]; xlabel = "x", ylabel = "particle spacing / h₀", title = "Lagrangian particles at t = 0.2 (N = 800, τ = 1e-5)")
ax2 = Axis(fig[1, 2]; xlabel = "x", ylabel = "ρ", title = "density on the same particles")
for (i, c) in ((1, :blue), (2, :black))
    d = studio_data(i); m, v = first(d.methods)
    x = vec(v.points); u = vec(v.u); p = sortperm(x); x, u = x[p], u[p]
    h0 = 1.0 / 800
    lines!(ax1, x[1:end-1] .+ diff(x) ./ 2, diff(x) ./ h0; color = c, label = string(m))
    lines!(ax2, x, u; color = c, label = string(m))
    println(rpad(string(m), 32), " spacing / h0: min ", round(minimum(diff(x)) / h0, digits = 2), ", max ", round(maximum(diff(x)) / h0, digits = 2))
end
xlims!(ax1, -0.5, 0.5); axislegend(ax1; position = :lt); vlines!(ax2, [0.185, 0.35]; color = :gray, linestyle = :dash)
save(joinpath(IMG, "e4c_particle_spacing.png"), fig)
println(joinpath(IMG, "e4c_particle_spacing.png"))
