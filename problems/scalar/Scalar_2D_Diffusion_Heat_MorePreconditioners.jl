using Penguin
using IterativeSolvers
using LinearSolve
using Statistics
using CSV
using DataFrames
using SparseArrays
using LinearAlgebra

# Preconditioner packages
using ILUZero
using IncompleteLU
using AlgebraicMultigrid
using Preconditioners

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

#–– 2) Define Solver Wrappers with Preconditioners ––

# 1. ILUZero (ILU(0))
function gmres_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# 2. IncompleteLU (ILU(τ))
function gmres_ilu_tau(A, b; kwargs...)
    # τ = 0.01 is a common starting point
    P = IncompleteLU.ilu(A, τ=0.01)
    return gmres(A, b; Pl=P, kwargs...)
end

# 3. AlgebraicMultigrid (AMG - Ruge-Stuben)
function gmres_amg_rs(A, b; kwargs...)
    ml = AlgebraicMultigrid.ruge_stuben(A)
    P = AlgebraicMultigrid.aspreconditioner(ml)
    return gmres(A, b; Pl=P, kwargs...)
end

# 4. AlgebraicMultigrid (AMG - Smoothed Aggregation)
function gmres_amg_sa(A, b; kwargs...)
    ml = AlgebraicMultigrid.smoothed_aggregation(A)
    P = AlgebraicMultigrid.aspreconditioner(ml)
    return gmres(A, b; Pl=P, kwargs...)
end

# 5. Preconditioners.jl - Diagonal
function gmres_diag(A, b; kwargs...)
    P = Preconditioners.DiagonalPreconditioner(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# 6. Preconditioners.jl - Cholesky (Incomplete Cholesky)
# Note: Requires A to be SPD.
function gmres_ichol(A, b; kwargs...)
    try
        P = Preconditioners.CholeskyPreconditioner(A, 2)
        return gmres(A, b; Pl=P, kwargs...)
    catch e
        @warn "CholeskyPreconditioner failed (matrix might not be SPD). Fallback to Identity." exception=e
        return gmres(A, b; kwargs...)
    end
end

methods = Dict(
  "GMRES (No Precond)" => gmres,
  "GMRES + ILU0"       => gmres_ilu0,
  "GMRES + ILU(0.01)"  => gmres_ilu_tau,
  "GMRES + AMG(RS)"    => gmres_amg_rs,
  "GMRES + AMG(SA)"    => gmres_amg_sa,
  "GMRES + Diagonal"   => gmres_diag,
  "GMRES + IChol(2)"   => gmres_ichol
)

#–– 3) Storage for residual histories ––
res_hist = Dict{String, Vector{Float64}}()

#–– 4) Loop over methods ––
for (name, meth) in methods
    println("→ Testing $name…")
    solver = DiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, u0, "BE")
    
    try
        solve_DiffusionUnsteadyMono!(
            solver, Fluide, Δt, Tend, bc_b, bc, "BE";
            method = meth,
            reltol   = eps(Float64),
            log    = true
        )
    catch err
        @warn "Solver failed for $name with atol=eps(). Retrying with abstol=1e-10..." exception=err
        solve_DiffusionUnsteadyMono!(
            solver, Fluide, Δt, Tend, bc_b, bc, "BE";
            method = meth,
            atol = 1e-10,
            btol = 1e-10,
            log    = true
        )
    end

    if !isempty(solver.ch)
        ch = solver.ch[1]
        resid = ch[:resnorm]
        niter = length(resid)
        final_res = resid[end]
        println("    iterations = $niter, final residual = $(final_res)")
        res_hist[name] = resid
    else
        println("    No convergence history found.")
        res_hist[name] = Float64[]
    end
end

#–– 5) Save to CSV ––
max_iter = isempty(res_hist) ? 0 : maximum(length, values(res_hist))
df = DataFrame(Iteration = 1:max_iter)

for (name, resid) in res_hist
    col_data = Vector{Union{Float64, Missing}}(missing, max_iter)
    col_data[1:length(resid)] = resid
    df[!, name] = col_data
end

results_dir = joinpath(BENCH_ROOT, "results", "scalar")
mkpath(results_dir)
output_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_MorePreconditioners.csv")
CSV.write(output_file, df)
println("Results saved to $output_file")
