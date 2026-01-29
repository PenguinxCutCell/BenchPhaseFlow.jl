# Scalar 2D Diffusion Heat AllPreconditioners Summary

No description available; please add a docstring to the originating problem script.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_AllPreconditioners_Summary.csv`

| Method | Iterations | FinalResidual | L2Error |
| --- | --- | --- | --- |
| GMRES (No Precond) | 118 | 1.01473e-16 | 0.00726905 |
| BICGSTAB + AMG(SA) | 17 | 3.78551e-16 | 0.00726905 |
| BICGSTAB + Diag+ILU0 | 25 | 1.6248e-14 | 0.00726905 |
| GMRES + AMG(RS) | 45 | 2.42497e-15 | 0.00726905 |
| GMRES + ILU0 | 11 | 2.75555e-15 | 0.00726905 |
| GMRES + Diagonal | 28 | 1.85232e-15 | 0.00726905 |
| GMRES + AMG(SA) | 58 | 1.48427e-15 | 0.00726905 |
| BICGSTAB + Diagonal | 8 | 6.29762e-17 | 0.00726905 |
| GMRES + Diag+ILU0 | 91 | 5.04356e-13 | 0.00726905 |
| GMRES + ILU(0.01) | 28 | 1.85232e-15 | 0.00726905 |
| BICGSTAB + ILU(0.01) | 8 | 1.28937e-16 | 0.00726905 |
| BICGSTAB + AMG(RS) | 20 | 2.84344e-16 | 0.00726905 |
| BICGSTAB (No Precond) | 27 | 2.64065e-17 | 0.00726905 |
| BICGSTAB + LLDL | 4 | 2.95441e-17 | 0.00726905 |
| GMRES + LLDL | 15 | 2.68927e-15 | 0.00726905 |
| BICGSTAB + ILU0 | 3 | 3.71342e-18 | 0.00726905 |
