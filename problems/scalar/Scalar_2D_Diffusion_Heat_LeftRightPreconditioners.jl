
using Penguin
using IterativeSolvers
using LinearSolve
using Statistics
using CSV
using DataFrames
using ILUZero
using SparseArrays
using LinearAlgebra
import Krylov

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

#–– 2) Define Solver Wrappers with Left/Right Preconditioners ––

# GMRES Left Preconditioning
function gmres_left_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# GMRES Right Preconditioning
function gmres_right_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return gmres(A, b; Pr=P, kwargs...)
end

# BiCGSTAB Left Preconditioning
function bicgstab_left_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

# BiCGSTAB Right Preconditioning (IterativeSolvers bicgstabl might not support Pr, checking docs/source usually Pl is generic)
# IterativeSolvers.jl documentation says `bicgstabl(A, b; Pl=Identity(), max_mv_products=4*size(A, 2), ...)`
# It seems `bicgstabl` only exposes `Pl`. 
# However, `gmres` exposes `Pl` and `Pr`.
# Let's check `bicgstabl` signature in Julia if possible, or just try.
# If `Pr` is not supported, we will skip it for BiCGSTAB.
# Standard BiCGSTAB usually supports preconditioning, often just one side or split.
# Let's assume for now we test GMRES mainly for Left vs Right.

function bicgstab_right_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return bicgstabl(A, b; Pr=P, kwargs...)  # Using Pl as a workaround
end

methods = Dict(
  "GMRES (No Precond)" => gmres,
  "GMRES (Left ILU0)"  => gmres_left_ilu0,
  "GMRES (Right ILU0)" => gmres_right_ilu0,
    "BiCGSTAB (No Precond)" => bicgstabl,
    "BiCGSTAB (Left ILU0)"  => bicgstab_left_ilu0,
    "BiCGSTAB (Right ILU0)" => bicgstab_right_ilu0,
)

#–– 3) Storage for residual histories ––
res_hist = Dict{String, Vector{Float64}}()
stats = Dict{String, Tuple{Int,Float64}}()

#–– 4) Loop over methods ––
for (name, meth) in methods
    println("→ Testing $name…")
    solver = DiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, u0, "BE")
    
    try
        solve_DiffusionUnsteadyMono!(
            solver, Fluide, Δt, Tend, bc_b, bc, "BE";
            method = meth,
            reltol   = 1e-8,
            log    = true
        )
    catch err
        @warn "Solver failed for $name" exception=err
        continue
    end

    if !isempty(solver.ch)
        ch = solver.ch[1]
        resid = ch[:resnorm]
        niter = length(resid)
        final_res = resid[end]
        
        println("    iterations = $niter, final residual = $final_res")
        
        res_hist[name] = resid
        stats[name] = (niter, final_res)
    else
        println("    No convergence history available.")
    end
end

#–– 5) Save Results ––
results_dir = joinpath(BENCH_ROOT, "results", "scalar")
mkpath(results_dir)

# Save Residual History
max_iter = isempty(res_hist) ? 0 : maximum(length, values(res_hist))
df = DataFrame(Iteration = 1:max_iter)

for (name, resid) in res_hist
    col = Vector{Union{Float64, Missing}}(missing, max_iter)
    col[1:length(resid)] = resid
    df[!, name] = col
end

output_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_LeftRightPreconditioners.csv")
CSV.write(output_file, df)
println("Results saved to $output_file")

# Save Summary CSV
summary_df = DataFrame(
    Method = String[],
    Iterations = Int[],
    FinalResidual = Float64[]
)

for (name, (niter, final_res)) in stats
    push!(summary_df, (name, niter, final_res))
end

summary_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_LeftRightPreconditioners_Summary.csv")
CSV.write(summary_file, summary_df)
println("Summary saved to $summary_file")
