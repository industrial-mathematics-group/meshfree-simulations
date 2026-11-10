# E5: particle management — the particle count at t = 0.2 relative to the initial N, over N, for the three variants with
# adaptivity (τ = 1e-5 and 1e-2), with the statistic final_particle_count on the U axis.
include(joinpath(@__DIR__, "session.jl"))
adapt = [m for m in ALL_METHODS if occursin("adapt", string(m))]
studio_set!(:methods => adapt)
studio_set!(:plot_style => :lines, :x_axis => :Grid_Ns__1, :u_axis => :final_particle_count, :compare_target => :PDE_tau,
            :compare_columns => 3, :plot_width => 470, :plot_height => 320)
studio_option!(:labels, :x_label, "N"); studio_option!(:labels, :y_label, "particles at t = 0.2")
studio_option!(:x_axis_1d, :log_scale, true); studio_option!(:y_axis_1d, :log_scale, true)
studio_press!(:update_plot)
println(shot("e5_particle_count"))
for (i, τ) in enumerate((1e-5, 1e-3, 1e-2)), (m, v) in studio_data(i).methods
    N = vec(v.x); n = vec(v.u)
    println("τ = $τ  ", rpad(string(m), 36), join(("N=$(Int(N[k])): $(round(n[k] / N[k], digits = 3))" for k in (1, 5, 9)), "  "))
end
