# JohansenColella P1 Neumann Convergence

Johansen-Colella Problem 1 (Neumann cut): constant-coefficient Poisson equation
outside a star-shaped boundary defined by
    Omega = {(r, theta) : r >= 0.25 + 0.05 cos(6 theta)}.
We solve Delta phi = 7 r^2 cos(3 theta) with Dirichlet data taken from the exact
solution phi(r, theta) = r^4 cos(3 theta) on the outer box, and Neumann data on
 the cut boundary from the analytical gradient.

**CSV source:** `results/scalar/JohansenColella_P1_Neumann_Convergence.csv`

| method | h | lp_norm | inside_cells | inside_cells_by_dim | all_err | full_err | cut_err | empty_err | pair_order_all | pair_order_full | pair_order_cut | trunc_all | trunc_full | trunc_cut | trunc_pair_order_all | trunc_pair_order_full | trunc_pair_order_cut | trunc_max_all | trunc_max_full | trunc_max_cut | trunc_pair_order_max_all | trunc_pair_order_max_full | trunc_pair_order_max_cut |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JohansenColella_P1_Neumann | 0.0625 | L^1 | 183 | [10, 10] | 0.0255608 | 0.0240627 | 0.00149811 | 0 | NaN | NaN | NaN | 0.0304974 | 0.0242112 | 0.061771 | NaN | NaN | NaN | 4.92801 | 0.105611 | 4.92801 | NaN | NaN | NaN |
| JohansenColella_P1_Neumann | 0.03125 | L^1 | 778 | [20, 20] | 0.0127479 | 0.0119124 | 0.000835523 | 0 | 1.00368 | 1.01434 | 0.842389 | 0.00924253 | 0.00852006 | 0.0181496 | 1.72233 | 1.50674 | 1.76699 | 0.175701 | 0.0481477 | 0.175701 | 4.80981 | 1.13322 | 4.80981 |
| JohansenColella_P1_Neumann | 0.015625 | L^1 | 3190 | [39, 39] | 0.00628555 | 0.00618691 | 9.8638e-05 | 0 | 1.02015 | 0.945168 | 3.08246 | 0.00370468 | 0.00309114 | 0.0129994 | 1.31894 | 1.46272 | 0.481496 | 0.0737805 | 0.0229776 | 0.0737805 | 1.25181 | 1.06724 | 1.25181 |
| JohansenColella_P1_Neumann | 0.0078125 | L^1 | 12938 | [77, 77] | 0.00315378 | 0.00313101 | 2.27734e-05 | 0 | 0.994954 | 0.98259 | 2.1148 | 0.00205703 | 0.00111725 | 0.0154472 | 0.848785 | 1.46819 | -0.248899 | 0.177295 | 0.011223 | 0.177295 | -1.26484 | 1.03377 | -1.26484 |
| JohansenColella_P1_Neumann | 0.00390625 | L^1 | 52072 | [154, 154] | 0.00158292 | 0.00157698 | 5.93499e-06 | 0 | 0.994498 | 0.989462 | 1.94003 | 0.00123856 | 0.000404193 | 0.0145697 | 0.731897 | 1.46683 | 0.084375 | 0.112009 | 0.00554608 | 0.112009 | 0.662533 | 1.01692 | 0.662533 |
| JohansenColella_P1_Neumann | 0.00195312 | L^1 | 208958 | [308, 308] | 0.000793311 | 0.00079181 | 1.50067e-06 | 0 | 0.996628 | 0.99394 | 1.98364 | 0.00092336 | 0.000152032 | 0.0159371 | 0.423701 | 1.41067 | -0.129422 | 0.743925 | 0.00593443 | 0.743925 | -2.73154 | -0.0976399 | -2.73154 |
