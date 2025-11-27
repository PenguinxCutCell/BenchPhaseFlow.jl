using Penguin
using IterativeSolvers
using LinearSolve
using Statistics
using CSV
using DataFrames
using ILUZero
using SparseArrays

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

# Wrapper for GMRES with ILU0
function gmres_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# Wrapper for BiCGSTAB with ILU0
function bicgstabl_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

# Wrapper for CG with ILU0 (Note: CG usually requires SPD and Cholesky, but we test if ILU helps or if it runs)
# ILU0 produces LU, not LLT. Preconditioned CG requires P to be SPD. 
# P = LU is not necessarily SPD. 
# We will focus on GMRES and BiCGSTAB which handle general matrices.

methods = Dict(
  "GMRES"         => gmres,
  "GMRES+ILU0"    => gmres_ilu0,
  "BiCGSTAB"      => bicgstabl,
  "BiCGSTAB+ILU0" => bicgstabl_ilu0,
  "CG"            => cg
)

#–– 3) Storage for residual histories ––
res_hist = Dict{String, Vector{Float64}}()

#–– 4) Loop over methods ––
for (name, meth) in methods
    println("→ Testing $name…")
    # Re-initialize solver for each method to ensure clean state
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
        # Retry with looser tolerance if strict tolerance fails
        solve_DiffusionUnsteadyMono!(
            solver, Fluide, Δt, Tend, bc_b, bc, "BE";
            method = meth,
            atol = 1e-10,
            btol = 1e-10,
            log    = true
        )
    end

    # Extract history from the first timestep (or the last one executed)
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
output_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_Preconditioners.csv")
CSV.write(output_file, df)
println("Results saved to $output_file")
