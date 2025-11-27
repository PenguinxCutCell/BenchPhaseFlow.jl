using Penguin
using IterativeSolvers
using LinearSolve
using Statistics
using CSV
using DataFrames
using SpecialFunctions
using Roots

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

#–– 1) Problem setup (identical to your Heat.jl) ––
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

#–– 2) The list of iterative‐solver methods to try ––
methods = Dict(
  "CG"      => cg,
  "BiCGSTAB"=> bicgstabl,
  "MINRES"  => minres,
  "IDRS"    => idrs,
  "GMRES"   => gmres,
   #"LSMR"    => lsmr,
  "LSQR"    => lsqr
)

#–– 3) Storage for residual histories (first time‐step only) ––
res_hist = Dict{String, Vector{Float64}}()
stats = Dict{String, Tuple{Int,Float64}}()  # (niter, final_res)
errors = Dict{String, Float64}()

#–– 4) Loop over methods ––
solvers = Dict{String, Any}()

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
            atol = eps(Float64),
            btol = eps(Float64),
            log    = true
        )
    end

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
    solvers[name] = solver
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
output_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_IterativeSolvers.csv")
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

summary_file = joinpath(results_dir, "Scalar_2D_Diffusion_Heat_IterativeSolvers_Summary.csv")
CSV.write(summary_file, summary_df)
println("Summary saved to $summary_file")
