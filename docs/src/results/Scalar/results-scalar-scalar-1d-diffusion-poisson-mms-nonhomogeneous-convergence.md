# Scalar 1D Diffusion Poisson MMS NonHomogeneous Convergence

1D MMS Poisson with non-homogeneous Dirichlet boundary conditions.
Analytical solution:
    u(x) = u_left + (u_right - u_left) * (x-dx)/(lx-dx) + A * sin(pi*(x-dx)/(lx-dx))

**CSV source:** `results/scalar/Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence | 0.2 | L^2 | 5 | [5] | 0.0167693 | 0.0167693 | 0 | 0 | NaN | NaN | NaN |
| Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence | 0.1 | L^2 | 10 | [10] | 0.00342658 | 0.00342658 | 0 | 0 | 2.29099 | 2.29099 | NaN |
| Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence | 0.05 | L^2 | 20 | [20] | 0.00078618 | 0.00078618 | 0 | 0 | 2.12384 | 2.12384 | NaN |
| Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence | 0.025 | L^2 | 40 | [40] | 0.000188837 | 0.000188837 | 0 | 0 | 2.05772 | 2.05772 | NaN |
| Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence | 0.0125 | L^2 | 80 | [80] | 4.63044e-05 | 4.63044e-05 | 0 | 0 | 2.02792 | 2.02792 | NaN |
| Scalar_1D_Diffusion_Poisson_MMS_NonHomogeneous_Convergence | 0.00625 | L^2 | 160 | [160] | 1.14664e-05 | 1.14664e-05 | 0 | 0 | 2.01374 | 2.01374 | NaN |
