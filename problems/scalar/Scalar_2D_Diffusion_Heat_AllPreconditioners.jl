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

#–– 2) Define Solver Wrappers with Preconditioners ––

# 1. ILUZero (ILU(0))
function gmres_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# 2. IncompleteLU (ILU(τ))
function gmres_ilu_tau(A, b; kwargs...)
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

# 6. LimitedLDLFactorizations (LLDL)
function gmres_lldl(A, b; kwargs...)
    # LLDL requires A to be symmetric (or Hermitian)
    # We check if it is approximately symmetric
    if !issymmetric(A)
        # Try to symmetrize or just warn?
        # For this benchmark, we'll try to use it anyway if it runs, 
        # or fallback if it errors.
        # LLDL might expect a symmetric matrix type or just check values.
    end
    
    try
        F = lldl(A)
        # Make SPD if needed (optional for GMRES but good for stability)
        F.D .= abs.(F.D) 
        return gmres(A, b; Pl=F, kwargs...)
    catch e
        @warn "LLDL failed (matrix might not be symmetric or suitable). Fallback to Identity." exception=e
        return gmres(A, b; kwargs...)
    end
end

# 7. LinearSolve.ComposePreconditioner (Diagonal + ILU0 - just for testing composition)
function gmres_composed(A, b; kwargs...)
    P1 = Preconditioners.DiagonalPreconditioner(A)
    P2 = ilu0(A)
    P = LinearSolve.ComposePreconditioner(P1, P2)
    return gmres(A, b; Pl=P, kwargs...)
end

# BiCGSTAB with same preconditioners
function bicgstabl_ilu0(A, b; kwargs...)
    P = ilu0(A)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

function bicgstabl_ilu_tau(A, b; kwargs...)
    P = IncompleteLU.ilu(A, τ=0.01)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

function bicgstabl_amg_rs(A, b; kwargs...)
    ml = AlgebraicMultigrid.ruge_stuben(A)
    P = AlgebraicMultigrid.aspreconditioner(ml)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

function bicgstabl_amg_sa(A, b; kwargs...)
    ml = AlgebraicMultigrid.smoothed_aggregation(A)
    P = AlgebraicMultigrid.aspreconditioner(ml)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

function bicgstabl_diag(A, b; kwargs...)
    P = Preconditioners.DiagonalPreconditioner(A)
    return bicgstabl(A, b; Pl=P, kwargs...)
end

function bicgstabl_lldl(A, b; kwargs...)
    try
        F = lldl(A)
        F.D .= abs.(F.D)
        return bicgstabl(A, b; Pl=F, kwargs...)
    catch e
        @warn "LLDL failed (matrix might not be symmetric or suitable). Fallback to Identity." exception=e
        return bicgstabl(A, b; kwargs...)
    end
end

function bicgstabl_composed(A, b; kwargs...)
    P1 = Preconditioners.DiagonalPreconditioner(A)
    P2 = ilu0(A)
    P = LinearSolve.ComposePreconditioner(P1, P2)
    return bicgstabl(A, b; Pl=P, kwargs...)
end


methods = Dict(
  "GMRES (No Precond)" => gmres,
  "GMRES + ILU0"       => gmres_ilu0,
  "GMRES + ILU(0.01)"  => gmres_ilu_tau,
  "GMRES + AMG(RS)"    => gmres_amg_rs,
  "GMRES + AMG(SA)"    => gmres_amg_sa,
  "GMRES + Diagonal"   => gmres_diag,
  "GMRES + LLDL"       => gmres_lldl,
  "GMRES + Diag+ILU0"  => gmres_composed,
  "BICGSTAB (No Precond)" => bicgstabl,
  "BICGSTAB + ILU0"       => bicgstabl_ilu0,
    "BICGSTAB + ILU(0.01)"  => bicgstabl_ilu_tau,
    "BICGSTAB + AMG(RS)"    => bicgstabl_amg_rs,
    "BICGSTAB + AMG(SA)"    => bicgstabl_amg_sa,
    "BICGSTAB + Diagonal"   => bicgstabl_diag,
    "BICGSTAB + LLDL"       => bicgstabl_lldl,
    "BICGSTAB + Diag+ILU0"  => bicgstabl_composed
)

#–– 3) Storage for residual histories ––
res_hist = Dict{String, Vector{Float64}}()
stats = Dict{String, Tuple{Int,Float64}}()
errors = Dict{String, Float64}()

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
        
        # Compute L2 error
        include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))
        _, _, global_err, _, _, _ = check_convergence(u_analytical, solver, capacity, 2)
        
        println("    iterations = $niter, final residual = $(final_res), L2 error = $global_err")
        res_hist[name] = resid
        stats[name] = (niter, final_res)
        errors[name] = global_err
    else
        println("    No convergence history found.")
        res_hist[name] = Float64[]
        stats[name] = (0, NaN)
        errors[name] = NaN
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
output_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_AllPreconditioners.csv")
CSV.write(output_file, df)
println("Results saved to $output_file")

# Save Summary CSV
summary_df = DataFrame(
    Method = String[],
    Iterations = Int[],
    FinalResidual = Float64[],
    L2Error = Float64[]
)

for (name, (niter, final_res)) in stats
    push!(summary_df, (name, niter, final_res, errors[name]))
end

summary_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_AllPreconditioners_Summary.csv")
CSV.write(summary_file, summary_df)
println("Summary saved to $summary_file")
