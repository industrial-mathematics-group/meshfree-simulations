# Run this script like: 
# julia --project=env_local startup.jl local
# julia --project=env_headless startup.jl headless

# 1. CLI argument (highest priority)
# 2. Environment variable JULIA_SIM_MODE
# 3. Fallback default
mode = if length(ARGS) > 0
    ARGS[1]
else
    get(ENV, "JULIA_SIM_MODE", "headless")
end


using Logging
using ThreadPinning
using HyperCloud
using StaticArrays

try
    pinthreads(:cores)
catch e
    @warn "Could not pin threads to cores!"
end

# Shared setup
min_level_to_show = Logging.Info
global_logger(ConsoleLogger(stderr, min_level_to_show))
PDEStudioCore.set_target_module!(Main)

if mode == "web"
    @eval using WGLMakie
    @eval using PDEStudio
    @eval using Bonito
    PDEStudioCore.set_save_path!("/bigtmp/pajung")
    HyperCloud.set_progress_interval!(0.5)
    HyperCloud.enable_progress_bar!(true)
    
elseif mode == "local"
    @eval using PDEStudio
    @eval using GLMakie
    PDEStudioCore.set_save_path!(pwd())
    HyperCloud.set_progress_interval!(0.5)
    HyperCloud.enable_progress_bar!(true)
    
elseif mode == "headless"
    PDEStudioCore.set_save_path!("/bigtmp/pajung")
    PDEStudioCore.set_stat_preset!("hyperbolic")
    HyperCloud.set_progress_interval!(60)
    HyperCloud.enable_progress_bar!(false)
else
    error("Unknown mode: $mode")
end