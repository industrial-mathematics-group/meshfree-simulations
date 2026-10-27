# Euler shock problems (data/review/euler/data)

### mach3step, 900x300 CFL 0.8 Float64  (Δmass relative to the final mass of RK2MUSCL2(EPD1,AR,MLS2))

| Method | t_end | min ρ | min p | mass change | Δmass | steps | runtime |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1(none,RU,MLS2) | 4 | 8.714e-01 | 1.000e+00 | +0.8028 | +0.2 % | 13166 | 616 s |
| RK2MUSCL2(EPD1,AR,MLS2) | 4 | 2.028e-01 | 1.860e-01 | +0.7990 | +0.0 % | 13582 | 928 s |
| RK2MUSCL2(EPD2,AR,MLS2) | 4 | 2.933e-01 | 2.742e-01 | +0.8039 | +0.3 % | 13528 | 1274 s |
| RK2MUSCL2(EPD2,RU,MLS2) | 4 | 2.944e-01 | 2.750e-01 | +0.8040 | +0.3 % | 13529 | 1253 s |
| RK2MUSCL2(VK,AR,MLS2) | 4 | 2.197e-01 | 2.333e-01 | +0.7989 | -0.0 % | 13533 | 776 s |
| RK3MUSCL3(EPD1,AR,MLS0) | **0.00339 (failed)** | -1.827e+17 | -4.427e+17 | -487794934737.6342 | -27115317515294.9 % | 33 | 6 s |
| RK3MUSCL3(EPD1,AR,MLS2) | 4 | 1.976e-01 | 1.933e-01 | +0.7975 | -0.1 % | 13598 | 1847 s |
| RK3MUSCL3(EPD1,AR,MLS3) | 4 | 2.312e-01 | 2.144e-01 | +0.7677 | -1.7 % | 40334 | 7815 s |
| RK3MUSCL3(EPD2,AR,MLS0) | **0.00358 (failed)** | -4.320e+18 | -1.429e+19 | -9233589971168.9004 | -513272496379167.7 % | 37 | 8 s |
| RK3MUSCL3(EPD2,AR,MLS2) | 4 | 2.895e-01 | 2.817e-01 | +0.8043 | +0.3 % | 13531 | 2488 s |
| RK3MUSCL3(EPD2,AR,MLS3) | 4 | 1.913e-01 | 1.680e-01 | +0.7634 | -2.0 % | 40559 | 9754 s |
| RK3MUSCL3(EPD2,RU,MLS2) | 4 | 2.888e-01 | 2.805e-01 | +0.8044 | +0.3 % | 13532 | 2340 s |
| RK4SSPMUSCL4(EPD1,AR,MLS0) | **0.00385 (failed)** | -2.659e+06 | -8.206e+06 | -0.6844 | -82.5 % | 38 | 14 s |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 4 | 6.588e-02 | 6.835e-02 | +0.7948 | -0.2 % | 13660 | 5117 s |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 4 | 1.377e-01 | 1.426e-01 | +0.7619 | -2.1 % | 41300 | 17909 s |
| RK4SSPMUSCL4(EPD2,AR,MLS0) | **0.00761 (failed)** | -1.165e+18 | -4.093e+18 | -2087237844087.1853 | -116024404605048.5 % | 77 | 30 s |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 4 | 2.577e-01 | 2.645e-01 | +0.8012 | +0.1 % | 13559 | 5140 s |
| RK4SSPMUSCL4(EPD2,RU,MLS2) | 4 | 2.481e-01 | 2.676e-01 | +0.8012 | +0.1 % | 13557 | 4960 s |

### rectangular, 1440x450 CFL 0.8 Float64  (Δmass relative to the final mass of RK2MUSCL2(EPD1,AR,MLS2))

| Method | t_end | min ρ | min p | mass change | Δmass | steps | runtime |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1(none,RU,MLS-) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0370 | -0.2 % | 3965 | 496 s |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0447 | +0.0 % | 4174 | 750 s |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0435 | -0.0 % | 4140 | 1085 s |
| RK2MUSCL2(EPD2,RU,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0435 | -0.0 % | 4139 | 1082 s |
| RK2MUSCL2(VK,AR,MLS2) | 0.2 | 1.398e+00 | 9.983e-01 | +3.0444 | -0.0 % | 4209 | 655 s |
| RK3MUSCL3(EPD1,AR,MLS0) | **0.0054 (failed)** | -6.961e+06 | -1.373e+08 | -0.2146 | -80.6 % | 304 | 138 s |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0450 | +0.0 % | 4224 | 1437 s |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.2 | 1.399e+00 | 9.762e-01 | +3.0418 | -0.1 % | 12187 | 6180 s |
| RK3MUSCL3(EPD2,AR,MLS0) | **0.00423 (failed)** | -2.415e+07 | -4.990e+08 | -8.2319 | -278.8 % | 242 | 134 s |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0438 | -0.0 % | 4139 | 2042 s |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.2 | 1.399e+00 | 1.000e+00 | +3.0441 | -0.0 % | 12177 | 8237 s |
| RK3MUSCL3(EPD2,RU,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0438 | -0.0 % | 4137 | 2087 s |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0454 | +0.0 % | 4292 | 4487 s |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0441 | -0.0 % | 12214 | 15143 s |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0437 | -0.0 % | 4134 | 4457 s |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.2 | 1.400e+00 | 9.894e-01 | +3.0443 | -0.0 % | 12184 | 16906 s |
| RK4SSPMUSCL4(EPD2,RU,MLS2) | 0.2 | 1.400e+00 | 1.000e+00 | +3.0437 | -0.0 % | 4144 | 4502 s |

### rectangular, 960x300 CFL 0.8 Float64  (Δmass relative to the final mass of RK2MUSCL2(EPD1,AR,MLS2))

| Method | t_end | min ρ | min p | mass change | Δmass | steps | runtime |
|---|---|---|---|---|---|---|---|
| RK4SSPMUSCL4(EPD1,AR,MLS0) | **0.00706 (failed)** | -4.990e+47 | -1.020e+49 | -30939082159377598336098240641710017216512.0000 | – | 265 | 120 s |
| RK4SSPMUSCL4(EPD2,AR,MLS0) | **0.145 (failed)** | -3.229e+19 | -6.000e+21 | +76529126064981.8594 | – | 5984 | 3454 s |

### rectangular/riemann, 800x800 CFL 0.8 Float64  (Δmass relative to the final mass of RK2MUSCL2(EPD1,AR,MLS2))

| Method | t_end | min ρ | min p | mass change | Δmass | steps | runtime |
|---|---|---|---|---|---|---|---|
| RK2MUSCL1(none,RU,MLS2) | 0.3 | 1.380e-01 | 2.900e-02 | +0.3848 | -0.0 % | 1644 | 193 s |
| RK2MUSCL2(EPD1,AR,MLS2) | 0.3 | 1.380e-01 | 2.597e-02 | +0.3853 | +0.0 % | 1815 | 318 s |
| RK2MUSCL2(EPD2,AR,MLS2) | 0.3 | 1.380e-01 | 2.877e-02 | +0.3850 | -0.0 % | 1758 | 458 s |
| RK2MUSCL2(EPD2,RU,MLS2) | 0.3 | 1.380e-01 | 2.889e-02 | +0.3850 | -0.0 % | 1758 | 447 s |
| RK2MUSCL2(VK,AR,MLS2) | 0.3 | 1.347e-01 | 2.701e-02 | +0.3854 | +0.0 % | 1874 | 276 s |
| RK3MUSCL3(EPD1,AR,MLS0) | **0.000824 (failed)** | -8.167e+14 | -4.491e+14 | -1326400415.8779 | -95746152641.2 % | 16 | 8 s |
| RK3MUSCL3(EPD1,AR,MLS2) | 0.3 | 1.379e-01 | 2.143e-02 | +0.3854 | +0.0 % | 1818 | 602 s |
| RK3MUSCL3(EPD1,AR,MLS3) | 0.3 | 1.377e-01 | 1.820e-02 | +0.3846 | -0.0 % | 5301 | 2788 s |
| RK3MUSCL3(EPD2,AR,MLS0) | **0.000983 (failed)** | -1.598e+05 | -7.662e+04 | -0.2040 | -42.5 % | 19 | 11 s |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.3 | 1.380e-01 | 2.896e-02 | +0.3850 | -0.0 % | 1758 | 872 s |
| RK3MUSCL3(EPD2,AR,MLS3) | 0.3 | 1.380e-01 | 2.526e-02 | +0.3852 | -0.0 % | 5253 | 3627 s |
| RK3MUSCL3(EPD2,RU,MLS2) | 0.3 | 1.380e-01 | 2.895e-02 | +0.3850 | -0.0 % | 1757 | 860 s |
| RK4SSPMUSCL4(EPD1,AR,MLS0) | **0.00073 (failed)** | -5.839e+15 | -3.229e+15 | -9443895253.8721 | -681707141716.1 % | 14 | 14 s |
| RK4SSPMUSCL4(EPD1,AR,MLS2) | 0.3 | 1.380e-01 | 2.743e-02 | +0.3855 | +0.0 % | 1834 | 1918 s |
| RK4SSPMUSCL4(EPD1,AR,MLS4) | 0.3 | 1.380e-01 | 2.201e-02 | +0.3844 | -0.1 % | 5324 | 6797 s |
| RK4SSPMUSCL4(EPD2,AR,MLS0) | **0.000988 (failed)** | -2.982e+04 | -1.652e+04 | -0.0149 | -28.9 % | 19 | 21 s |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.3 | 1.380e-01 | 2.899e-02 | +0.3850 | -0.0 % | 1769 | 1821 s |
| RK4SSPMUSCL4(EPD2,AR,MLS4) | 0.3 | 1.380e-01 | 2.784e-02 | +0.3857 | +0.0 % | 5272 | 7474 s |
| RK4SSPMUSCL4(EPD2,RU,MLS2) | 0.3 | 1.380e-01 | 2.900e-02 | +0.3850 | -0.0 % | 1767 | 1846 s |

### rectangular/riemann/bc=Dict(3 => :outflow, 2 => :outflow, 1 => :outflow, 4 => :outflow, 5 => :outflow)-outflow, 800x800 CFL 0.8 Float64  (Δmass relative to the final mass of RK2MUSCL2(EPD1,AR,MLS2))

| Method | t_end | min ρ | min p | mass change | Δmass | steps | runtime |
|---|---|---|---|---|---|---|---|
| RK2MUSCL2(EPD2,AR,MLS2) | 0.3 | 1.380e-01 | 2.897e-02 | +0.3788 | – | 1760 | 494 s |
| RK3MUSCL3(EPD2,AR,MLS2) | 0.3 | 1.380e-01 | 2.894e-02 | +0.3787 | – | 1758 | 954 s |
| RK4SSPMUSCL4(EPD2,AR,MLS2) | 0.3 | 1.380e-01 | 2.900e-02 | +0.3732 | – | 1768 | 1997 s |
