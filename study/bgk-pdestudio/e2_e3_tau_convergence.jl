# E2: the relaxation time τ (1e-5, 1e-3, 1e-2) at N = 800 for three variants.
# E3: convergence of the L1 density error at t = 0.2 over N, per τ, all nine variants (statistic on the U axis).
include(joinpath(@__DIR__, "session.jl"))
three = [Symbol("ARS222MUSCL2(VK) static"), Symbol("ARS222MUSCL2(EPD2) lagrangian"), Symbol("ARS233MUSCL3(EPD2) lagrangian")]
studio_set!(:methods => three)
studio_set!(:Grid_Ns__1 => 800, :t => 0.2, :compare_target => :PDE_tau, :compare_columns => 3, :plot_style => :particles_lines,
            :plot_width => 470, :plot_height => 300)
studio_option!(:labels, :x_label, "x"); studio_option!(:labels, :y_label, "ρ"); studio_option!(:plot_style, :line_width, 1.2)
studio_press!(:update_plot)
println(shot("e2_tau_N800"))

studio_set!(:methods => ALL_METHODS)
studio_set!(:plot_style => :lines, :x_axis => :Grid_Ns__1, :u_axis => :l1error, :t => 0.2, :compare_target => :PDE_tau,
            :compare_columns => 3, :plot_width => 470, :plot_height => 340)
studio_option!(:labels, :x_label, "N"); studio_option!(:labels, :y_label, "L1 error of ρ")
studio_option!(:x_axis_1d, :log_scale, true); studio_option!(:y_axis_1d, :log_scale, true)
studio_press!(:update_plot)
println(studio_state().status.text)
println(shot("e3_l1_convergence"))
# observed orders between N = 400 and 1600 per variant and τ
for (i, τ) in enumerate((1e-5, 1e-3, 1e-2))
    d = studio_data(i)
    for (m, v) in d.methods
        N = vec(v.x); e = vec(v.u)
        i4, i16 = argmin(abs.(N .- 400)), argmin(abs.(N .- 1600))
        println("τ = $τ  ", rpad(string(m), 36), " L1(1600) = ", round(e[i16], sigdigits = 3), "  order 400→1600 = ", round(log(e[i4] / e[i16]) / log(4), digits = 2))
    end
end
