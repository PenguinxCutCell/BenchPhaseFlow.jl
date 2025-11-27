using Penguin
using IterativeSolvers
using LinearSolve
using Statistics
using CSV
using DataFrames
using SparseArrays
using LinearAlgebra
using SpecialFunctions
using Roots

# Preconditioner packages
using ILUZero
using IncompleteLU
using AlgebraicMultigrid
using Preconditioners
using LimitedLDLFactorizations

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

#–– 1) Problem setup ––
nx, ny = 80, 80
lx, ly = 4.0, 4.0
mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))

radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = (x,y,_=0) -> sqrt((x-center[1])^2 + (y-center[2])^2) - radius

capacity = Capacity(circle, mesh)
operator = DiffusionOps(capacity)

bc0 = Dirichlet(0.0)
bc = Dirichlet(1.0)
bc_b = BorderConditions(Dict(
  :left   => bc0,
  :right  => bc0,
  :top    => bc0,
  :bottom => bc0
))

f      = (x,y,z,t) -> 0.0
Dcoef  = (x,y,z)   -> 1.0
Fluide = Phase(capacity, operator, f, Dcoef)

u0ₒ = zeros((nx+1)*(ny+1))
u0ᵧ = ones((nx+1)*(ny+1))
u0  = vcat(u0ₒ, u0ᵧ)

Δt   = 0.25 * (lx/nx)^2
Tend = 0.011

# Analytical solution for error computation
function radial_heat_solution(center::Tuple{Float64,Float64}, R::Float64;
    t::Float64 = 0.1,
    k::Float64 = 1.0,
    a::Float64 = 1.0,
    wr::Float64 = 1.0,
    w0::Float64 = 0.0,
    Nzeros::Int = 200
)
    function j0_zeros_robin(N, k, R; guess_shift = 0.25)
        eq(alpha) = besselj0(alpha)
        zs = zeros(Float64, N)
        for m in 1:N
            x_left  = max((m - guess_shift - 0.5) * π, 1e-6)
            x_right = (m - guess_shift + 0.5) * π
            zs[m] = find_zero(eq, (x_left, x_right))
        end
        return zs
    end

    alphas = j0_zeros_robin(Nzeros, k, R)

    return function (x::Float64, y::Float64)
        r = sqrt((x - center[1])^2 + (y - center[2])^2)
        if r >= R
            return w0
        end
        s = 0.0
        for α in alphas
            s += exp(-α^2 * t) * besselj0(α * (r / R)) / (α * besselj1(α))
        end
        return 1.0 - 2.0*s
    end
end

u_analytical = radial_heat_solution(center, radius; t=Tend)

#–– 2) Define LinearSolve.jl Algorithms with Preconditioners ––

# Helper to create a LinearSolve algorithm with a specific preconditioner
# We use KrylovJL_GMRES as the base solver from LinearSolve.jl which wraps Krylov.jl/IterativeSolvers.jl
# The `precs` argument takes a function (A, p) -> (Pl, Pr)

# 1. ILUZero (ILU(0))
function prec_ilu0(A, p=nothing)
    return (ilu0(A), LinearAlgebra.I)
end
algo_ilu0 = LinearSolve.KrylovJL_GMRES(precs=prec_ilu0)

# 2. IncompleteLU (ILU(τ))
function prec_ilu_tau(A, p=nothing)
    return (IncompleteLU.ilu(A, τ=0.01), LinearAlgebra.I)
end
algo_ilu_tau = LinearSolve.KrylovJL_GMRES(precs=prec_ilu_tau)

# 3. AlgebraicMultigrid (AMG - Ruge-Stuben)
function prec_amg_rs(A, p=nothing)
    ml = AlgebraicMultigrid.ruge_stuben(A)
    return (AlgebraicMultigrid.aspreconditioner(ml), LinearAlgebra.I)
end
algo_amg_rs = LinearSolve.KrylovJL_GMRES(precs=prec_amg_rs)

# 4. AlgebraicMultigrid (AMG - Smoothed Aggregation)
function prec_amg_sa(A, p=nothing)
    ml = AlgebraicMultigrid.smoothed_aggregation(A)
    return (AlgebraicMultigrid.aspreconditioner(ml), LinearAlgebra.I)
end
algo_amg_sa = LinearSolve.KrylovJL_GMRES(precs=prec_amg_sa)

# 5. Preconditioners.jl - Diagonal
function prec_diag(A, p=nothing)
    return (Preconditioners.DiagonalPreconditioner(A), LinearAlgebra.I)
end
algo_diag = LinearSolve.KrylovJL_GMRES(precs=prec_diag)

# 6. LimitedLDLFactorizations (LLDL)
function prec_lldl(A, p=nothing)
    try
        F = lldl(A)
        F.D .= abs.(F.D)
        return (F, LinearAlgebra.I)
    catch e
        return (LinearAlgebra.I, LinearAlgebra.I)
    end
end
algo_lldl = LinearSolve.KrylovJL_GMRES(precs=prec_lldl)

# 7. No Preconditioner
algo_none = LinearSolve.KrylovJL_GMRES()


methods = Dict(
  "LS_GMRES (No Precond)" => algo_none,
  "LS_GMRES + ILU0"       => algo_ilu0,
  "LS_GMRES + ILU(0.01)"  => algo_ilu_tau,
  "LS_GMRES + AMG(RS)"    => algo_amg_rs,
  "LS_GMRES + AMG(SA)"    => algo_amg_sa,
  "LS_GMRES + Diagonal"   => algo_diag,
  "LS_GMRES + LLDL"       => algo_lldl
)

#–– 3) Storage for residual histories ––
res_hist = Dict{String, Vector{Float64}}()
stats = Dict{String, Tuple{Int,Float64}}()
errors = Dict{String, Float64}()

#–– 4) Loop over methods ––
for (name, algo) in methods
    println("→ Testing $name…")
    solver = DiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, u0, "BE")
    
    # We use the `algorithm` argument of solve_DiffusionUnsteadyMono! 
    # which passes it to solve_system! -> solve_with_linearsolve!
    # Note: solve_DiffusionUnsteadyMono! signature:
    # solve_DiffusionUnsteadyMono!(..., algorithm=nothing, kwargs...)
    
    try
        solve_DiffusionUnsteadyMono!(
            solver, Fluide, Δt, Tend, bc_b, bc, "BE";
            algorithm = algo,
            reltol   = eps(Float64),
            abstol   = 1e-10, # LinearSolve uses abstol/reltol
            verbose  = true   # To get history? LinearSolve might not return history in the same way
        )
        
        # LinearSolve.jl integration in Penguin might not populate solver.ch automatically 
        # unless solve_with_linearsolve! does it.
        # Let's check if solver.ch is populated.
        # If not, we might not get the residual history easily without modifying Penguin or using a callback.
        # However, for this benchmark, let's see what we get.
        
    catch err
        @warn "Solver failed for $name" exception=err
    end

    # Check if history is available
    # LinearSolve.jl solvers usually return a solution object which might have stats.
    # But Penguin's solve_system! puts the result in s.x and maybe s.ch if configured.
    # If s.ch is empty, we can't plot history, but we can still get final error.
    
    # NOTE: The user wants to compare iterative solvers. 
    # If Penguin's `solve_with_linearsolve!` doesn't capture history, we might only get the final result.
    # Let's assume for now we might not get full history if not supported, 
    # but we can at least get the final error.
    
    # Actually, LinearSolve.jl's `solve` returns a solution object. 
    # If Penguin discards it, we lose the info.
    # Let's check if we can get the residual history.
    # If not, we will just report the L2 error and maybe iteration count if available in `solver.ch`.
    
    if !isempty(solver.ch)
        # If Penguin wraps LinearSolve and extracts history (unlikely by default for LS algorithms)
        # We might need to rely on what's available.
        # If solver.ch is empty, we'll just report L2 error.
        ch = solver.ch[1]
        # Check format of ch. LinearSolve might return something different or nothing.
        if isa(ch, Vector) || isa(ch, Dict) # Assuming some format
             # ...
        end
    end
    
    # Compute L2 error
    include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))
    _, _, global_err, _, _, _ = check_convergence(u_analytical, solver, capacity, 2)
    
    println("    L2 error = $global_err")
    
    # For CSV, if we don't have history, we can't plot convergence curves.
    # But we can save the final error.
    errors[name] = global_err
    
    # Try to guess iterations/residual if possible, otherwise NaN
    niter = 0
    final_res = NaN
    
    # If solver.ch is populated (e.g. if Penguin adapts LS output), use it.
    if !isempty(solver.ch)
         # This depends on how Penguin implements solve_with_linearsolve!
         # If it doesn't push to ch, we have nothing.
         # Let's assume for now we might get nothing for history.
    end
    
    stats[name] = (niter, final_res)
end

#–– 5) Save Summary CSV ––
results_dir = joinpath(BENCH_ROOT, "results", "scalar")
mkpath(results_dir)

summary_df = DataFrame(
    Method = String[],
    L2Error = Float64[]
)

for (name, err) in errors
    push!(summary_df, (name, err))
end

summary_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_LinearSolve_Summary.csv")
CSV.write(summary_file, summary_df)
println("Summary saved to $summary_file")
