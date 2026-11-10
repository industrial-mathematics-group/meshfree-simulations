# E1: the Sod structure of all nine variants (3 methods × static / Lagrangian / Lagrangian + particle management)
# at t = 0.2, N = 800, τ = 1e-5 (near the Euler limit), and a zoom on the contact and the shock.
include(joinpath(@__DIR__, "session.jl"))
studio_set!(:Grid_Ns__1 => 800, :PDE_tau => 1e-5, :t => 0.2, :compare_target => :methods, :compare_columns => 3,
            :plot_style => :particles_lines, :plot_width => 470, :plot_height => 260)
studio_option!(:labels, :x_label, "x"); studio_option!(:labels, :y_label, "ρ")
studio_option!(:plot_style, :line_width, 1.5)
studio_press!(:update_plot)
println(shot("e1_methods_N800_tau1e-5"))
studio_view!(xlims = (0.05, 0.45), ylims = (0.1, 0.5))
println(shot("e1_methods_N800_tau1e-5_zoom"))
# the plateau values between rarefaction and contact (x ≈ 0.05) and between contact and shock (x ≈ 0.27)
for i in 1:9
    d = studio_data(i); m, v = first(d.methods)
    x = vec(v.points); u = vec(v.u)
    at(x0) = u[argmin(abs.(x .- x0))]
    println(rpad(string(m), 36), " ρ(0.05) = ", round(at(0.05), digits = 4), "  ρ(0.27) = ", round(at(0.27), digits = 4))
end
