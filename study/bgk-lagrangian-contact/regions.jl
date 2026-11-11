# L1 density error of 1-D Sod runs split into regions (against the exact Euler solution at t = 0.2), the shock position
# (where ρ crosses (ρ*_R + ρ_R)/2, exact 0.3504) and the change of the mass Σ h_i ρ_i (h_i: particle spacing) from t = 0.
# Regions: rarefaction [-0.30, -0.01], left plateau [-0.01, 0.16], contact [0.16, 0.21], right plateau [0.21, 0.33],
# shock [0.33, 0.37], rest.
#     julia --project=env_web study/bgk-lagrangian-contact/regions.jl CONFIG SAVE_DIR N METHOD...
include(joinpath(@__DIR__, "common.jl"))
global_logger(ConsoleLogger(stderr, Logging.Error))
name, save, N = ARGS[1], ARGS[2], parse(Int, ARGS[3])

const REGIONS = (("rarefaction", -0.30, -0.01), ("left plateau", -0.01, 0.16), ("contact", 0.16, 0.21),
                 ("right plateau", 0.21, 0.33), ("shock", 0.33, 0.37))
@printf("exact: contact %.4f  shock %.4f  ρ*L %.5f  ρ*R %.5f\n", 0.2ust, 0.2S, ρsL, ρsR)
println("| Variant | ", join((r[1] for r in REGIONS), " | "), " | rest | total | shock x | mass change |")
println("|---|", repeat("---|", length(REGIONS) + 4))
for m in ARGS[4:end]
    d = load_run(name, save, m, N); s = snapshot(d); s0 = snapshot(d, 1)
    e = abs.(s.ρ .- ρexact.(s.x, d.t[end])) .* s.h
    parts = [sum(e[i] for i in eachindex(s.x) if a <= s.x[i] < b; init = 0.0) for (_, a, b) in REGIONS]
    total = sum(e)
    dm = sum(s.h .* s.ρ) / sum(s0.h .* s0.ρ) - 1
    @printf("| %s | %s | %.2e | %.2e | %.4f | %+.2e |\n", m, join((@sprintf("%.2e", v) for v in parts), " | "),
            total - sum(parts), total, shock_position(s), dm)
end
