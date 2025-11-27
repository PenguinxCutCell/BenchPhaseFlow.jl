
using Penguin
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using CSV
using DataFrames
using KrylovPreconditioners
using RandomizedPreconditioners
import Krylov
using SpecialFunctions
using Roots

const BENCH_ROOT = joinpath(@__DIR__, "..", "..")

#–– 1) Define the Problem ––
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

# Analytical solution (for error check)
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


#–– 2) Define Preconditioner Wrappers ––

# Wrapper to handle SubArray inputs for Preconditioners that don't support them
struct KPWrapper{P}
    precond::P
end

function LinearAlgebra.ldiv!(y::AbstractVector, W::KPWrapper, x::AbstractVector)
    x_vec = Vector(x)
    y_vec = similar(x_vec)
    ldiv!(y_vec, W.precond, x_vec)
    copy!(y, y_vec)
    return y
end

function LinearAlgebra.ldiv!(W::KPWrapper, x::AbstractVector)
    x_vec = Vector(x)
    ldiv!(W.precond, x_vec)
    copy!(x, x_vec)
    return x
end

# KrylovPreconditioners.jl
# ILU0
function prec_kp_ilu0(A)
    # ilu(A) returns ILU0Preconditioner which supports ldiv!
    P = KrylovPreconditioners.ilu(A, τ=0.0) 
    return KPWrapper(P)
end

# Wrappers for IterativeSolvers
function gmres_kp_ilu0(A, b; kwargs...)
    P = prec_kp_ilu0(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# Krylov.jl Solver Wrappers
function krylov_gmres(A, b; kwargs...)
    atol = get(kwargs, :abstol, 1e-8)
    rtol = get(kwargs, :reltol, 1e-8)
    x, stats = Krylov.gmres(A, b; atol=atol, rtol=rtol, history=true)
    return x, stats.residuals
end

function krylov_bicgstab(A, b; kwargs...)
    atol = get(kwargs, :abstol, 1e-8)
    rtol = get(kwargs, :reltol, 1e-8)
    x, stats = Krylov.bicgstab(A, b; atol=atol, rtol=rtol, history=true)
    return x, stats.residuals
end

function krylov_cg(A, b; kwargs...)
    atol = get(kwargs, :abstol, 1e-8)
    rtol = get(kwargs, :reltol, 1e-8)
    x, stats = Krylov.cg(A, b; atol=atol, rtol=rtol, history=true)
    return x, stats.residuals
end

methods = Dict(
  "GMRES + KP_ILU0"        => gmres_kp_ilu0,
  "Krylov.gmres"           => krylov_gmres,
  "Krylov.bicgstab"        => krylov_bicgstab,
  "Krylov.cg"              => krylov_cg
)

# Wrappers for IterativeSolvers
function gmres_kp_ilu0(A, b; kwargs...)
    P = prec_kp_ilu0(A)
    return gmres(A, b; Pl=P, kwargs...)
end

# Krylov.jl Solver Wrapper
# Krylov.jl solvers return (x, stats)
function krylov_gmres(A, b; kwargs...)
    # Map kwargs to Krylov.jl options if needed
    # Krylov.jl uses 'atol', 'rtol', 'itmax', 'history'
    # IterativeSolvers uses 'abstol', 'reltol', 'maxiter', 'log'
    
    # Extract tolerances
    atol = get(kwargs, :abstol, 1e-8)
    rtol = get(kwargs, :reltol, 1e-8)
    
    x, stats = Krylov.gmres(A, b; atol=atol, rtol=rtol, history=true)
    
    # Return in IterativeSolvers format: (x, history)
    # We need to construct a history object or just return the residual vector
    # Penguin expects (x, history) where history is a vector of residuals or similar
    # Actually Penguin checks `solver.ch` which is populated by the method return if it's standard.
    # But here we are passing the method to `solve_DiffusionUnsteadyMono!`.
    # Penguin's `solve_linear_system` handles the return.
    # If we return (x, stats), Penguin might not know how to parse stats.
    # Let's look at how Penguin handles it.
    # Penguin usually expects `IterativeSolvers` return type if `log=true`.
    
    # We will return (x, stats.residuals) to mimic IterativeSolvers history vector
    return x, stats.residuals
end

methods = Dict(
  "GMRES + KP_ILU0"        => gmres_kp_ilu0,
  "Krylov.gmres"           => krylov_gmres,
  "Krylov.bicgstab"        => krylov_bicgstab,
  "Krylov.cg"              => krylov_cg
)

#–– 3) Run Benchmark ––
res_hist = Dict{String, Vector{Float64}}()
stats = Dict{String, Tuple{Int,Float64}}()
errors = Dict{String, Float64}()

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
        # Penguin stores the second return value of the solver method in solver.ch
        # For IterativeSolvers, it's a ConvergenceHistory object.
        # For our Krylov wrapper, it's a Vector.
        
        history_obj = solver.ch[1]
        
        resid = Float64[]
        if isa(history_obj, Vector)
            resid = history_obj
        elseif hasproperty(history_obj, :data) && haskey(history_obj.data, :resnorm)
             resid = history_obj[:resnorm]
        else
             # Fallback
             println("    Unknown history format: $(typeof(history_obj))")
        end
        
        if !isempty(resid)
            niter = length(resid)
            final_res = resid[end]
            
            # Compute L2 error
            include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))
            _, _, global_err, _, _, _ = check_convergence(u_analytical, solver, capacity, 2)
            
            println("    iterations = $niter, final residual = $final_res, L2 error = $global_err")
            
            res_hist[name] = resid
            stats[name] = (niter, final_res)
            errors[name] = global_err
        end
    else
        println("    No convergence history available.")
    end
end

#–– 4) Save Results ––
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

output_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_KrylovPreconditioners.csv")
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

summary_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_KrylovPreconditioners_Summary.csv")
CSV.write(summary_file, summary_df)
println("Summary saved to $summary_file")
