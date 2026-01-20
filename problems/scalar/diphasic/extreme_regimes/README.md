# Extremes regimes for the diphasic heat equation

This folder contains benchmarks for extreme regimes of the diphasic heat equation, where the diffusivity ratio and henry jump are very large or very small.

## 1D Diphasic Heat Equation

$\frac{\partial T^\pm}{\partial t} = D^\pm \frac{\partial^2 T^\pm}{\partial x^2}$
with interface conditions at $x=x_\Gamma$:
- $T^+ = He \cdot T^-$
- $D^+ \frac{\partial T^+}{\partial x} = D^- \frac{\partial T^-}{\partial x}$
### Parameters
- Domain: $x \in [0, L]$ with $L=20.0$
- Interface position: $x_\Gamma = 10.01$
- Time interval: $t \in [0, T_{end}]$ with $T_{end} = 0.1$
- Initial condition: $T(x,0) = 0.0$
- Boundary conditions: $T(0,t) = 1.0$, $T(L,t) = 0.0$
- Diffusivity in phase -: $D^- = 1.0$
- Diffusivity in phase +: $D^+ = D_{ratio} \cdot D^-$ with $D_{ratio} \in \{1.0, 100.0, 10000.0, 1.0^6\}$
- Henry jump: $He \in \{1.0, 10.0, 100.0\}$ 

| case                             |   all |  full |   cut |
| :------------------------------- | ----: | ----: | ----: | 
| Heat_2ph_1D_ratio1.0_He1.0       | 1.655 | 1.499 | 2.695 |  
| Heat_2ph_1D_ratio1.0_He10.0      | 1.655 | 1.499 | 2.671 | 
| Heat_2ph_1D_ratio1.0_He100.0     | 1.655 | 1.499 | 2.665 | 
| Heat_2ph_1D_ratio100.0_He1.0     | 1.632 | 1.477 | 2.806 |  
| Heat_2ph_1D_ratio100.0_He10.0    | 1.623 | 1.473 | 2.849 |  
| Heat_2ph_1D_ratio100.0_He100.0   | 1.622 | 1.472 | 2.855 | 
| Heat_2ph_1D_ratio10000.0_He1.0   | 1.154 | 1.005 | 1.898 |
| Heat_2ph_1D_ratio10000.0_He10.0  | 1.151 | 1.002 | 1.929 |
| Heat_2ph_1D_ratio10000.0_He100.0 | 1.151 | 1.002 | 1.929 |
| Heat_2ph_1D_ratio1.0e6_He1.0     | 1.606 | 1.453 | 2.436 |
| Heat_2ph_1D_ratio1.0e6_He10.0    | 1.606 | 1.457 | 2.457 |
| Heat_2ph_1D_ratio1.0e6_He100.0   | 1.606 | 1.457 | 2.459 |  
*Table: Convergence rates for various extreme regimes of the diphasic heat equation in 1D. 'all' refers to the overall convergence rate, 'full' to the convergence rate in the full cells, and 'cut' to the convergence rate in the cut cells near the interface.*

| case | ‖e‖_all | ‖e‖_full | ‖e‖_cut | phase-1 (all) | phase-2 (all) |
|---|---:|---:|---:|---:|---:|
| Heat_2ph_1D_ratio1.0_He1.0 | 4.56e-4 | 4.56e-4 | 1.2e-5 | 4.56e-4 | 4.56e-4 |
| Heat_2ph_1D_ratio1.0_He10.0 | 8.27e-4 | 8.27e-4 | 2.3e-5 | 8.27e-4 | 8.27e-4 |
| Heat_2ph_1D_ratio1.0_He100.0 | 9.00e-4 | 9.00e-4 | 2.5e-5 | 9.00e-4 | 9.00e-4 |
| Heat_2ph_1D_ratio1.0e6_He1.0 | 9.93e-4 | 9.91e-4 | 6.9e-5 | 9.80e-4 | 9.93e-4 |
| Heat_2ph_1D_ratio1.0e6_He10.0 | 9.81e-4 | 9.78e-4 | 6.9e-5 | 9.81e-4 | 9.17e-4 |
| Heat_2ph_1D_ratio1.0e6_He100.0 | 9.81e-4 | 9.78e-4 | 6.9e-5 | 9.81e-4 | 9.10e-4 |
| Heat_2ph_1D_ratio100.0_He1.0 | 8.60e-4 | 8.60e-4 | 1.2e-5 | 6.80e-4 | 8.60e-4 |
| Heat_2ph_1D_ratio100.0_He10.0 | 9.04e-4 | 9.04e-4 | 1.0e-5 | 7.40e-4 | 9.04e-4 |
| Heat_2ph_1D_ratio100.0_He100.0 | 9.09e-4 | 9.09e-4 | 1.0e-5 | 7.46e-4 | 9.09e-4 |
| Heat_2ph_1D_ratio10000.0_He1.0 | 8.13e-3 | 8.11e-3 | 5.7e-4 | 8.13e-3 | 1.73e-3 |
| Heat_2ph_1D_ratio10000.0_He10.0 | 8.21e-3 | 8.19e-3 | 5.7e-4 | 8.21e-3 | 9.68e-4 |
| Heat_2ph_1D_ratio10000.0_He100.0 | 8.22e-3 | 8.20e-3 | 5.8e-4 | 8.22e-3 | 9.15e-4 |
*Table: Errors for various extreme regimes of the diphasic heat equation in 1D at the finest resolution (nx=256). 'all' refers to the overall error, 'full' to the error in the full cells, and 'cut' to the error in the cut cells near the interface. The last two columns show the errors in phase-1 and phase-2 respectively.*

Observations:   
- The overall convergence rates remain close to 1.5-1.6 for moderate diffusivity ratios (1.0 and 100.0) but drop to around 1.15 for very high diffusivity ratios (10000.0 and 1.0e6).
- Errors remains relatively low even in extreme regimes, indicating the robustness of the cut-cell numerical method.

## 2D Diphasic Heat Equation (extreme regimes)

- `Heat_2ph_2D.jl`: sweeps Henry jumps `[0.01, 0.1, 1.0, 10.0, 100.0]` at the default diffusivity ratio.
- `Heat_2ph_2D_ratio.jl`: sweeps diffusivity ratios `[1.0, 100.0, 10000.0, 1.0e6]` at a fixed Henry jump (default `He=2.0`, configurable via `He` argument).

Each benchmark writes one CSV per case under `results/scalar/diphasic/extreme_regimes/` by default.

# 2D Diphasic Heat Equation
$\frac{\partial T^\pm}{\partial t} = \nabla \cdot (D^\pm \nabla T^\pm)$
with interface conditions at the circular interface
- $T^+ = He \cdot T^-$
- $D^+ \nabla T^+ \cdot n = D^- \nabla T^- \cdot n$

### Parameters

- Domain: $(x,y) \in [0, L_x] \times [0, L_y]$
- Interface: Circle centered at $(x_c, y_c)$ with radius $R$
- Time interval: $t \in [0, T_{end}]$ with $T_{end} = 0.1$
- Initial condition: $T(x,y,0) = 0.0$
- Boundary conditions: $T(x=0,y,t) = 1.0$, $T(x=L_x,y,t) = 0.0$
- Diffusivity in phase -: $D^- = 1.0$
- Diffusivity in phase +: $D^+ = D_{ratio} \cdot D^-$ with $D_{ratio} \in \{1.0, 100.0, 10000.0, 1.0^6\}$
- Henry jump: $He \in \{0.01, 0.1, 1.0, 10.0, 100.0\}$

|   He | p_fit(all_err) | p_fit(full_err) | p_fit(cut_err) | err_full | err_cut |
| ---: | -------------: | --------------: | -------------: | -------------------------------: | ------------------------------: |
| 0.01 |          1.830 |           1.651 |          2.442 |                         2.14e-05 |                        1.35e-06 |
|  0.1 |          1.790 |           1.616 |          2.416 |                         1.92e-04 |                        1.14e-05 |
|    1 |          1.642 |           1.496 |          2.310 |                         9.70e-04 |                        5.50e-05 |
|   10 |          1.544 |           1.422 |          2.088 |                         1.66e-03 |                        2.11e-04 |
|  100 |          1.479 |           1.408 |          1.972 |                         2.00e-03 |                        2.94e-04 |
*Table: Convergence rates and errors for various Henry jumps (He) in the 2D diphasic heat equation with diffusivity ratio D_ratio=1.0 at the finest resolution (nx=128, ny=128). 'all_err', 'full_err', and 'cut_err' refer to the overall, full cell, and cut cell errors respectively.*

Observations:
- Convergence rates remain robust across a wide range of Henry jumps, with overall rates ranging from approximately 1.48 to 1.83.
- Errors still very low even for extreme Henry jumps, demonstrating the effectiveness of the numerical method in handling sharp interface conditions.
- Working even for a few cells inside the circle.

| ratio | p_all | p_full | p_cut | h_finest | all_err_finest | cut_err_finest | N_int
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.607 | 1.453 | 2.294 | 0.0625 | 3.23e-3 | 2.95e-4 | 32
| 100 | 1.539 | 1.371 | 2.001 | 0.0625 | 3.74e-3 | 6.85e-4 | 32
| 1e4 | 1.540 | 1.369 | 1.993 | 0.0625 | 3.72e-3 | 7.06e-4 | 32
| 1e6 | 1.475 | 1.301 | 1.952 | 0.0625 | 4.89e-3 | 8.30e-4 | 32
*Table: Convergence rates and errors for various diffusivity ratios (D_ratio) in the 2D diphasic heat equation with Henry jump He=2.0 at the finest resolution (nx=128, ny=128). 'all_err_finest' and 'cut_err_finest' refer to the overall and cut cell errors respectively. 32 cells per diameter.*

Observations:
- Convergence rates slightly decrease with increasing diffusivity ratio, particularly in the full cells.
- Errors remain low even for extreme diffusivity ratios, indicating the robustness of the numerical method in handling large contrasts in material properties.