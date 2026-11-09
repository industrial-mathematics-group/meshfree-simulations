#!/bin/bash
# The three options against U2 as on main: (1) limited fallback MOODFallbackLimiter(none → VK/BJ), (2) VK + MOODPAD,
# (3) EPD1 instead of EPD2. CPU, Float64, reduced sizes.
cd "$(dirname "$0")/../.."   # repo root; run with an env that develops the HyperCloudCore branch for the fallback variants (ENV_PROJECT)
r() { local tag=$1; shift; echo -n "[$tag] "; timeout 60m julia -t 8 --project=${ENV_PROJECT:-env_headless} tools/mood_u2_shear/e2e_u2.jl "$@" T=Float64 2>&1 | grep -E "steps=|ERROR|argmin" | tail -2 | tr '\n' ' '; echo; }
variants() {   # CASE TMAX [grid]
  local c=$1 t=$2; shift 2
  r K3-EPD2       $c "RK3MUSCL3(EPD2,AR,MLS2)" $t "$@"
  r K3-EPD1       $c "RK3MUSCL3(EPD1,AR,MLS2)" $t "$@"
  r K3-EPD2+VKfb  $c "RK3MUSCL3(EPD2,AR,MLS2)" $t "$@" Limiter_fallback=:VK
  r K3-EPD1+VKfb  $c "RK3MUSCL3(EPD1,AR,MLS2)" $t "$@" Limiter_fallback=:VK
  r K3-EPD2+BJfb  $c "RK3MUSCL3(EPD2,AR,MLS2)" $t "$@" Limiter_fallback=:BJ
  r K4-EPD2       $c "RK4SSPMUSCL4(EPD2,AR,MLS2)" $t "$@"
  r K4-EPD2+VKfb  $c "RK4SSPMUSCL4(EPD2,AR,MLS2)" $t "$@" Limiter_fallback=:VK
  r K2-EPD2       $c "RK2MUSCL2(EPD2,AR,MLS2)" $t "$@"
  r K2-EPD1       $c "RK2MUSCL2(EPD1,AR,MLS2)" $t "$@"
  r K2-VK         $c "RK2MUSCL2(VK,AR,MLS2)" $t "$@"
  r K2-VK+PAD     $c "RK2MUSCL2(VK,AR,MLS2)" $t "$@" MOOD_criterion=:PAD
  r K3-VK+PAD     $c "RK3MUSCL3(EPD1,AR,MLS2)" $t "$@" MOOD_criterion=:PAD Limiter_name=:VK Limiter_mode=:hard
}
variants KelvinHelmholtz 0.5 "Grid_Ns=(128,128)"
variants Gresho 1.0 "Grid_Ns=(100,100)"
variants IsentropicVortex_Float64 10.0 "Grid_Ns=(100,100)"
variants Sod1D 0.2
variants Riemann2D_C4 0.25 "Grid_Ns=(150,150)"
variants Riemann2D_C12 0.25 "Grid_Ns=(150,150)"
variants Implosion 1.0 "Grid_Ns=(100,100)"
variants Noh_Float64 0.3 "Grid_Ns=(100,100)"
echo DONE
