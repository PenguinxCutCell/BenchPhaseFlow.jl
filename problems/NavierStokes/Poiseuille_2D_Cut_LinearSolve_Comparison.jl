using Pkg
project_path = joinpath(@__DIR__, "..", "..")
Pkg.activate(project_path)

# Check and add packages
required_pkgs = ["LinearSolve", "Krylov", "IterativeSolvers", "ILUZero"]
installed_pkgs = Pkg.project().dependencies
for pkg in required_pkgs
    if !haskey(installed_pkgs, pkg)
        println("Installing $pkg...")
        Pkg.add(pkg)
    end
end

using Penguin
using SparseArrays
using LinearAlgebra
using LinearSolve
using Krylov
using IterativeSolvers
using ILUZero
using Statistics
using CSV
using DataFrames
using Test

# --- Global Stats Storage ---
const LINEAR_SOLVER_STATS = Dict{String, Any}()

function reset_stats!()
    empty!(LINEAR_SOLVER_STATS)
    LINEAR_SOLVER_STATS["total_iters"] = 0
    LINEAR_SOLVER_STATS["total_time"] = 0.0
    LINEAR_SOLVER_STATS["calls"] = 0
    LINEAR_SOLVER_STATS["history"] = []
end

# --- Monkey Patch Penguin.solve_navierstokes_linear_system! ---
# We need to use Penguin internals.
function Penguin.solve_navierstokes_linear_system!(s::Penguin.NavierStokesMono; method=Base.:\, algorithm=nothing, kwargs...)
    # Access internal function via Penguin.
    Ared, bred, keep_idx_rows, keep_idx_cols = Penguin.remove_zero_rows_cols!(s.A, s.b)

    kwargs_nt = (; kwargs...)
    precond_builder = haskey(kwargs_nt, :precond_builder) ? kwargs_nt.precond_builder : nothing
    if precond_builder !== nothing
        kwargs_nt = Base.structdiff(kwargs_nt, (precond_builder=precond_builder,))
    end

    precond_kwargs = (;)
    if precond_builder !== nothing
        precond_result = try
            precond_builder(Ared, s)
        catch err
            if err isa MethodError
                precond_builder(Ared)
            else
                rethrow(err)
            end
        end
        # We assume Penguin._preconditioner_kwargs exists and works, or we handle it manually if needed.
        # Since we verified it exists, we call it.
        precond_kwargs = Penguin._preconditioner_kwargs(precond_result)
    end

    solve_kwargs = merge(kwargs_nt, precond_kwargs)

    xred = nothing
    
    # Timing and Stats
    t0 = time()
    iters = 0
    
    if algorithm !== nothing
        prob = LinearSolve.LinearProblem(Ared, bred)
        # Force log=true to get stats if possible, though LinearSolve might not always return it in sol.iters
        # We rely on sol.iters if available.
        sol = LinearSolve.solve(prob, algorithm; solve_kwargs...)
        xred = sol.u
        
        # Try to extract iterations
        if hasproperty(sol, :iters)
            iters = sol.iters
        elseif hasproperty(sol, :stats) && hasproperty(sol.stats, :niter)
             iters = sol.stats.niter
        end
        
    elseif method === Base.:\
        try
            xred = Ared \ bred
            iters = 1 # Direct solve
        catch e
            if e isa SingularException
                @warn "Direct solver hit SingularException; falling back to bicgstabl" sizeA=size(Ared)
                xred = IterativeSolvers.bicgstabl(Ared, bred)
                iters = -1 # Unknown
            else
                rethrow(e)
            end
        end
    else
        log = get(solve_kwargs, :log, false)
        if log
            xred, ch = method(Ared, bred; solve_kwargs...)
            push!(s.ch, ch)
            iters = ch.iters
        else
            xred = method(Ared, bred; solve_kwargs...)
            iters = -1
        end
    end
    
    t1 = time()
    dt = t1 - t0
    
    # Update global stats
    if haskey(LINEAR_SOLVER_STATS, "total_iters")
        LINEAR_SOLVER_STATS["total_iters"] += iters
        LINEAR_SOLVER_STATS["total_time"] += dt
        LINEAR_SOLVER_STATS["calls"] += 1
        push!(LINEAR_SOLVER_STATS["history"], (iters, dt))
    end

    N = size(s.A, 2)
    s.x = zeros(N)
    s.x[keep_idx_cols] = xred
    return s
end

# --- Problem Setup ---
const nx, ny = 64, 64
const Lx, Ly = 2.0, 1.0
const x0, y0 = 0.0, 0.0

const y_wall_bot = 0.2
const y_wall_top = 0.8
const channel_height = y_wall_top - y_wall_bot

const Umax = 1.0
const μ = 1.0
const ρ = 1.0

body = (x, y, _=0) -> begin
    if y < y_wall_bot
        return y_wall_bot - y
    elseif y > y_wall_top
        return y - y_wall_top
    else
        return -min(y - y_wall_bot, y_wall_top - y)
    end
end

mesh_p = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
dx, dy = Lx / nx, Ly / ny
mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

capacity_ux = Capacity(body, mesh_ux; compute_centroids=false)
capacity_uy = Capacity(body, mesh_uy; compute_centroids=false)
capacity_p  = Capacity(body, mesh_p; compute_centroids=false)
operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

parabola = (x, y) -> begin
    if y < y_wall_bot || y > y_wall_top
        return 0.0
    else
        y_local = y - y_wall_bot
        return 4 * Umax * y_local * (channel_height - y_local) / (channel_height^2)
    end
end

ux_left  = Dirichlet(parabola)
ux_right = Dirichlet(parabola)
ux_bot   = Dirichlet((x, y)->0.0)
ux_top   = Dirichlet((x, y)->0.0)
bc_ux = BorderConditions(Dict(
    :left=>ux_left, :right=>ux_right, :bottom=>ux_bot, :top=>ux_top
))

uy_zero = Dirichlet((x, y)->0.0)
bc_uy = BorderConditions(Dict(:left=>uy_zero, :right=>uy_zero, :bottom=>uy_zero, :top=>uy_zero))

pressure_gauge = PinPressureGauge()
u_bc = Dirichlet(0.0)

fᵤ = (x, y, z=0.0) -> 0.0
fₚ = (x, y, z=0.0) -> 0.0

fluid = Fluid((mesh_ux, mesh_uy),
              (capacity_ux, capacity_uy),
              (operator_ux, operator_uy),
              mesh_p,
              capacity_p,
              operator_p,
              μ, ρ, fᵤ, fₚ)

nu = prod(operator_ux.size)
np = prod(operator_p.size)

# --- Preconditioner ---
function ilu_builder(A)
    return ILUZero.ilu0(A)
end

# --- Runner ---
function run_solver_test(solver_name, algo; use_precond=false)
    println("\nRunning solver: $solver_name")
    reset_stats!()
    
    x0_vec = zeros(4 * nu + np)
    solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, u_bc; x0=x0_vec)
    
    kwargs = Dict{Symbol, Any}(
        :tol => 1e-10,
        :maxiter => 50, # Reduced maxiter for safety
        :relaxation => 1.0,
        :algorithm => algo,
        :abstol => 1e-8,
        :reltol => 1e-8
    )
    
    if use_precond
        # We pass Pl directly to LinearSolve via kwargs if possible, 
        # OR we use Penguin's precond_builder if we want to use its mechanism.
        # But Penguin's mechanism calls _preconditioner_kwargs which might expect specific types.
        # Let's try passing Pl directly to LinearSolve via kwargs, 
        # but Penguin's solve_navierstokes_linear_system! merges kwargs.
        # So we can pass Pl = ...
        # But A changes every iteration? No, for steady Stokes it might be constant, but for Navier-Stokes it changes (Oseen).
        # So we need a builder.
        kwargs[:precond_builder] = ilu_builder
    end

    t_start = time()
    _, iters, res = solve_NavierStokesMono_steady!(solver; kwargs...)
    t_end = time()
    elapsed = t_end - t_start
    
    # Stats from linear solver
    total_linear_iters = LINEAR_SOLVER_STATS["total_iters"]
    avg_linear_iters = LINEAR_SOLVER_STATS["calls"] > 0 ? total_linear_iters / LINEAR_SOLVER_STATS["calls"] : 0
    total_linear_time = LINEAR_SOLVER_STATS["total_time"]

    println("  Time: $elapsed s")
    println("  Picard Iterations: $iters")
    println("  Residual: $res")
    println("  Total Linear Iters: $total_linear_iters")
    println("  Avg Linear Iters: $avg_linear_iters")

    # Compute errors
    uωx = solver.x[1:nu]
    xs = mesh_ux.nodes[1]
    ys = mesh_ux.nodes[2] .+ 0.5 * dy
    LIux = LinearIndices((length(xs), length(ys)))
    icol = Int(cld(length(xs), 2))
    ux_profile = [uωx[LIux[icol, j]] for j in 1:length(ys)]
    ux_analytical = [parabola(0.0, y) for y in ys]
    fluid_indices = findall(y -> y_wall_bot <= y <= y_wall_top, ys)
    profile_err = ux_profile[fluid_indices] .- ux_analytical[fluid_indices]
    ℓ2_profile = sqrt(sum(abs2, profile_err) / length(profile_err))
    ℓinf_profile = maximum(abs, profile_err[2:end-1])

    return (
        Solver = solver_name,
        Precond = use_precond,
        Time = elapsed,
        PicardIters = iters,
        Residual = res,
        TotalLinearIters = total_linear_iters,
        AvgLinearIters = avg_linear_iters,
        L2_Profile = ℓ2_profile,
        Linf_Profile = ℓinf_profile
    )
end

solvers = [
    ("UMFPACK", UMFPACKFactorization(), false),
    ("KrylovJL_GMRES+ILU0", KrylovJL_GMRES(), true),
    ("KrylovJL_BICGSTAB+ILU0", KrylovJL_BICGSTAB(), true),
    # ("IterativeSolversJL_GMRES", IterativeSolversJL_GMRES(), true), # Might need different precond setup
]

results = []
for (name, algo, use_precond) in solvers
    try
        res = run_solver_test(name, algo; use_precond=use_precond)
        push!(results, res)
    catch e
        println("Solver $name failed: $e")
        # push!(results, (Solver=name, Precond=use_precond, Time=NaN, PicardIters=-1, Residual=NaN, TotalLinearIters=-1, AvgLinearIters=-1, L2_Profile=NaN, Linf_Profile=NaN))
        showerror(stdout, e, catch_backtrace())
    end
end

df = DataFrame(results)
println("\nResults Summary:")
println(df)

csv_path = joinpath(@__DIR__, "..", "..", "results", "NavierStokes", "Poiseuille_2D_Cut_LinearSolve_Comparison.csv")
mkpath(dirname(csv_path))
CSV.write(csv_path, df)
println("Results written to $csv_path")
