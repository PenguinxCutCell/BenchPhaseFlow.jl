# A Space-Time Cartesian Cut-Cell Method for Two-Phase Diffusion Problems using a Two-Fluid Approach

Refer to the original publication for details:
[https://hal.science/hal-05457114](https://hal.science/hal-05457114)


## Benchmarks

Static problems :
Monophasic : 
- Johansen & Colella Problem 1 : 2D Star shaped poisson problem with Dirichlet BCs
- Johansen & Colella Problem 2 : 2D Flower-shaped laplace problem
- 2D Poisson equation in a disk with Robin BCs
- 3D Heat equation in a sphere with Robin BCs

Diphasic :
- Diphasic 2D unsteady diffusion with continuity/jump conditions at the interface
- Diphasic 2D unsteady diffusion with continuity/jump conditions at the interface with various Henry/diffusivity ratios
- Diphasic 3D unsteady diffusion with continuity conditions at the interface

Moving problems :
Monophasic :
- 2D Diffusion with an oscillating circular boundary
- 2D Diffusion with multiples ellipses (McCorquodale-Colella problem)
- 3D Schwartz-Colella diffusion on a expanding sphere

Diphasic :
- Manufactured solution with an oscillating two-phase interface

## Revision History

- 20/01/2026 : Paper send for review. Data and benchmarks using Penguin.jl v20/01/2026.
- 29/01/2026 : Updated benchmarks using Penguin.jl v29/01/2026 : Corrections in small cut-cell handling/time integration => Less oscillations in convergence results. Errors remains sensibly unchanged. Overall fit remain sensibly unchanged.

## Files

- problems/scalar/johansenColella/Problem1_PoissonConstant.jl : Johansen & Colella Problem 1 : 2D Star shaped poisson problem with Dirichlet BCs implementation
- problems/scalar/johansenColella/Problem2_FlowerLaplace.jl : Johansen & Colella Problem 2 : 2D Flower-shaped laplace problem implementation
- problems/scalar/Scalar_2D_Diffusion_Poisson_Robin.jl : 2D Poisson equation in a disk with Robin BCs implementation
- problems/scalar/Scalar_3D_Diffusion_Heat_Robin.jl : 3D Heat equation in a sphere with Robin BCs implementation
- problems/scalar/diphasic/Heat_2ph_2D.jl : Diphasic 2D unsteady diffusion with continuity/jump conditions at the interface implementation
- problems/scalar/diphasic/Hextreme_regimes/Heat_2ph_2D_ratio.jl, Heat_2ph_2D.jl, Heat_2ph_2D_ratio_he_dict.jl : Diphasic 2D unsteady diffusion with continuity/jump conditions at the interface with various Henry/diffusivity ratios implementation
- problems/scalar/diphasic/Heat_2ph_3D.jl : Diphasic 3D unsteady diffusion with continuity conditions at the interface implementation
- problems/scalar/PrescibedMotion/Heat_2D_Moving.jl, Heat_2D_Moving_CFL.jl : 2D Diffusion with an oscillating circular boundary implementation
- problems/scalar/PrescibedMotion/JohansenColella/MovingDirichlet.jl : 2D Diffusion with multiples ellipses (McCorquodale-Colella problem) implementation
- problems/scalar/PrescibedMotion/SchwartzColella3D/ExpandingSphere.jl : 3D Schwartz-Colella diffusion on a expanding sphere implementation
- problems/scalar/diphasic/PrescibedMotion/Heat_2ph_2D_MovingVertical.jl : Manufactured solution with an oscillating two-phase interface implementation

# Updated results 

Must do :
- heat 3d robin ok
- Diphasic 2d ok 
- diphasic 3d
- oscillating circle ok
- multiples ellipses ok
- expand sphere : in progress
- diphasic vertical oscillating ok

- Heat 3D Robin : Table 5 :
method,h,lp_norm,inside_cells,inside_cells_by_dim,all_err,full_err,cut_err,empty_err,pair_order_all,pair_order_full,pair_order_cut
Scalar_3D_Diffusion_Heat_Robin,2.0,L^2,0,[1, 1, 1],NaN,NaN,NaN,NaN,NaN,NaN,NaN
Scalar_3D_Diffusion_Heat_Robin,1.0,L^2,1,[2, 2, 2],0.3510253624198517,0.2962108088053846,0.18835594444786746,0.0,NaN,NaN,NaN
Scalar_3D_Diffusion_Heat_Robin,0.5,L^2,7,[4, 4, 4],0.1930974163822428,0.09073362508893155,0.1704524024228419,0.0,0.8622464111043169,1.7069150807447548,0.14409263142630127
Scalar_3D_Diffusion_Heat_Robin,0.25,L^2,147,[8, 8, 8],0.0584895706868391,0.0433281479429949,0.03928996659403885,0.0,1.7230775564942975,1.0663327299433407,2.1171360872589307
Scalar_3D_Diffusion_Heat_Robin,0.125,L^2,1599,[16, 16, 16],0.015406072229590527,0.0132598886752049,0.007843622496304452,0.0,1.9246803053902586,1.708235906881164,2.3245689339439592
Scalar_3D_Diffusion_Heat_Robin,0.0625,L^2,14915,[32, 32, 32],0.003919853643445592,0.0035991449440153428,0.0015528709727477396,0.0,1.9746274001590038,1.8813425547435172,2.3365821401868883
Scalar_3D_Diffusion_Heat_Robin,0.041666666666666664,L^2,52587,[48, 48, 48],0.0017484077968773844,0.001635596254950067,0.0006178630228253844,0.0,1.991167076209917,1.9451461201774176,2.2729303501270586
,,,,,,,,,,,
fit (all,full,cut) : (1.9030, 1.6465, 2.2756)

- Heat 2ph 2D Henry=1, ratio 1 : Table 6 : 

method,h,lp_norm,inside_cells,inside_cells_phase1,inside_cells_phase2,inside_cells_by_dim,all_err,full_err,cut_err,empty_err,phase1_all_err,phase1_full_err,phase1_cut_err,phase1_empty_err,phase2_all_err,phase2_full_err,phase2_cut_err,phase2_empty_err,pair_order_all,pair_order_full,pair_order_cut
Heat_2ph_2D_He1.0_Dl1.0,2.0,L^2,8,1,7,[1, 7],0.42756354912054156,0.2464342030431003,0.3494005897348569,0.0,0.42756354912054156,0.2464342030431003,0.3494005897348569,0.0,0.039966560899202315,0.01233886714116183,0.03801418613862994,0.0,NaN,NaN,NaN
Heat_2ph_2D_He1.0_Dl1.0,1.0,L^2,48,5,43,[5, 43],0.2084731229294724,0.14095124908791562,0.15360269647543065,0.0,0.2084731229294724,0.14095124908791562,0.15360269647543065,0.0,0.053059947771197,0.04431566286480379,0.02918013163670109,0.0,1.03627746641886,0.8060062412477856,1.1856785006998682
Heat_2ph_2D_He1.0_Dl1.0,0.5,L^2,224,37,187,[37, 187],0.07031683838012141,0.06442826394403635,0.028168361062257452,0.0,0.07031683838012141,0.06442826394403635,0.028168361062257452,0.0,0.023571278574386147,0.021336209932781406,0.010018548763947583,0.0,1.56791928856606,1.1294306367397595,2.4470560125854597
Heat_2ph_2D_He1.0_Dl1.0,0.25,L^2,960,177,783,[177, 783],0.0189039501744812,0.01857491044757759,0.0035116995975938347,0.0,0.0189039501744812,0.01857491044757759,0.0035116995975938347,0.0,0.006910983780701462,0.006730080487962969,0.0015708957453181962,0.0,1.8951824728943403,1.794338465718317,3.0038342825522215
Heat_2ph_2D_He1.0_Dl1.0,0.125,L^2,3968,749,3219,[749, 3219],0.004948875703371288,0.0048806451916200404,0.0008189464212845773,0.0,0.004948875703371288,0.0048806451916200404,0.0008189464212845773,0.0,0.0018670906738571956,0.0018312678202232076,0.00036398592695247874,0.0,1.933515019039665,1.9282114743115382,2.100328463994796
Heat_2ph_2D_He1.0_Dl1.0,0.0625,L^2,16128,3101,13027,[3101, 13027],0.0012320516992155057,0.0012215528951717378,0.00016049895276094568,0.0,0.0012320516992155057,0.0012215528951717378,0.00016049895276094568,0.0,0.0004650844256689466,0.0004611468963429266,6.039091813426955e-5,0.0,2.006038012226159,1.9983555399737711,2.351205184192987
fit (all,full,cut) : (1.9438, 1.9091, 2.4466)

- Heat 2ph 2D Henry : extreme regimes : DONE CF CSV

- Oscillating moving circle : Table 11 :

method,h,lp_norm,inside_cells,inside_cells_by_dim,all_err,full_err,cut_err,empty_err,pair_order_all,pair_order_full,pair_order_cut
Scalar_2D_Diffusion_Heat_Moving,1.0,L^2,5,[4, 4],1.232676848884725,1.0360684800145272,0.6678729808105446,0.0,NaN,NaN,NaN
Scalar_2D_Diffusion_Heat_Moving,0.5,L^2,21,[7, 7],0.6413763868665183,0.5398399487630786,0.34631849409095444,0.0,0.9425514959773699,0.940515715669635,0.9474743179808274
Scalar_2D_Diffusion_Heat_Moving,0.25,L^2,97,[13, 13],0.2492253322127605,0.09543930861768446,0.23022728897144062,0.0,1.3637205231626415,2.4998762448110585,0.5890405852854532
Scalar_2D_Diffusion_Heat_Moving,0.125,L^2,449,[26, 26],0.042743581999314405,0.03797507439371915,0.019619060296650524,0.0,2.5436709978759904,1.3295308008957916,3.5527310001864754
Scalar_2D_Diffusion_Heat_Moving,0.0625,L^2,1925,[51, 51],0.014166211099583187,0.013270378254190888,0.004957680698537497,0.0,1.593253867786749,1.516843297624715,1.9845186796073133
Scalar_2D_Diffusion_Heat_Moving,0.03125,L^2,7909,[102, 102],0.005373822149176953,0.005320090050253638,0.0007580279336390103,0.0,1.3984334659233193,1.3186869226104048,2.7093424386596
fit (all,full,cut) : (1.7935, 1.6176, 2.3209)

- Mutliples ellipses moving : Table 13 : 
method,h,lp_norm,inside_cells,inside_cells_by_dim,all_err,full_err,cut_err,empty_err,pair_order_all,pair_order_full,pair_order_cut
JohansenColella_Moving_Dirichlet,0.5,L^2,0,[1, 1],NaN,NaN,NaN,NaN,NaN,NaN,NaN
JohansenColella_Moving_Dirichlet,0.2222222222222222,L^2,0,[2, 2],0.002114456512118754,0.0,0.002114456512118754,0.0,NaN,NaN,NaN
JohansenColella_Moving_Dirichlet,0.125,L^2,2,[4, 4],0.0017645756788149203,0.00039261715618771115,0.0017203427260090388,0.0,0.31438795749456916,NaN,0.3585108365798596
JohansenColella_Moving_Dirichlet,0.06060606060606061,L^2,33,[7, 7],0.00024180037237862538,0.00020122930090927998,0.00013406785050117104,0.0,2.745546801596432,0.9232940009078924,3.525164220297699
JohansenColella_Moving_Dirichlet,0.031746031746031744,L^2,151,[12, 12],7.7184711357751e-5,6.261489832376488e-5,4.5131520861632036e-5,0.0,1.7659497499278323,1.8054326813455626,1.6837602331606962
JohansenColella_Moving_Dirichlet,0.015748031496062992,L^2,719,[24, 24],2.1746454275057748e-5,1.6616225947480735e-5,1.4028874110185941e-5,0.0,1.8069273391805336,1.8923316332423143,1.6667275331265081
JohansenColella_Moving_Dirichlet,0.0078125,L^2,3131,[49, 49],6.5323712457044825e-6,4.32204612271124e-6,4.898141627685286e-6,0.0,1.7156879552379376,1.9210688700909122,1.5011073356304039
,,,,,,,,,,,

- Diphasic 2D moving vertical interface : Table 15, 16, 17 : Done cf csv prescribed_motion/...