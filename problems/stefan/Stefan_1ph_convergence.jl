using Penguin
using LinearAlgebra
using SparseArrays
using SpecialFunctions
using Roots
using CSV
using DataFrames

# Mesh and time-step convergence study for one-phase Stefan problem.
# Varies mesh resolution and initial time, recording errors to CSV.

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

# ============================================================================
# Analytical solution functions
# ============================================================================

function find_lambda(Stefan_number)
    f = lambda -> sqrt(pi) * lambda * exp(lambda^2) * erf(lambda) - 1 / Stefan_number
    return find_zero(f, (1e-6, 10.0), Bisection())
end

function analytical_temperature(x, t, T0, k, lambda)
    t <= 0 && return T0
    return T0 - (T0 / erf(lambda)) * (erf(x / (2 * sqrt(k * t))))
end

analytical_position(t, k, lambda) = 2 * lambda * sqrt(k * t)

# ============================================================================
# Convergence study
# ============================================================================

# Test parameters
mesh_sizes = [8, 16, 32, 64, 128]
t_starts = [0.001, 0.01, 0.1, 0.2, 0.4]
delta_t_sim = 0.1  # Tend = Tstart + 0.1

# Stefan parameters (fixed)
Stefan_number = 1.0
lambda = find_lambda(Stefan_number)
lx = 1.0
x0 = 0.0

# Storage for results (per Tstart for CSVs)
results_by_tstart = Dict{Float64, Vector{NamedTuple}}()

println("Starting convergence study...")
println("Mesh sizes: ", mesh_sizes)
println("Starting times: ", t_starts, "\n")

# Loop over all combinations
total_runs = length(mesh_sizes) * length(t_starts)
run_count = 0

for nx in mesh_sizes
    for Tstart in t_starts
        global run_count += 1
        Tend = Tstart + delta_t_sim
        dx = lx / nx
        dt = 0.5 * dx^2

        println("[$run_count/$total_runs] nx=$nx, Tstart=$Tstart, Tend=$Tend, dx=$dx, dt=$dt")

        try
            # Build mesh and space-time mesh
            mesh = Penguin.Mesh((nx,), (lx,), (x0,))
            x_offset = mesh.nodes[1][1] - x0
            STmesh = Penguin.SpaceTimeMesh(mesh, [Tstart, Tstart + dt], tag=mesh.tag)

            # Initial interface position
            xf0_phys = analytical_position(Tstart, 1.0, lambda)
            xf = xf0_phys + x_offset
            body = (x, t, _=0) -> (x - xf)

            # Capacity and operators
            capacity = Capacity(body, STmesh)
            operator = DiffusionOps(capacity)

            # Boundary conditions
            bc = Dirichlet(0.0)
            bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(
                :top => Dirichlet(0.0),
                :bottom => Dirichlet(1.0)
            ))
            rho, L = 1.0, 1.0
            stef_cond = InterfaceConditions(nothing, FluxJump(1.0, 1.0, rho * L))

            # Source and diffusion coefficient
            f = (x, y, z, t) -> 0.0
            K = (x, y, z) -> 1.0

            # Phase definition
            fluid = Phase(capacity, operator, f, K)

            # Initial condition
            x_nodes = mesh.nodes[1]
            x_nodes_phys = x_nodes .- x_offset
            u0o = analytical_temperature.(x_nodes_phys, Tstart, 1.0, 1.0, lambda)
            u0y = zeros(nx + 1)
            u0 = vcat(u0o, u0y)

            # Solver
            solver = MovingLiquidDiffusionUnsteadyMono(fluid, bc_b, bc, dt, u0, mesh, "BE")

            # Solve using simplified direct method
            solver, xf_log, stefan_residuals = solve_MovingLiquidDiffusionUnsteadyMono_Simple!(
                solver, fluid, xf, dt, Tstart, Tend, bc_b, bc, stef_cond, mesh, "BE";
                method=Base.:\, max_inner_iter=1, tol=1e-8, damping=1.0
            )

            # Post-processing
            xf_log_phys = xf_log .- x_offset
            xf_num = xf_log_phys[end]
            xf_exact_Tend = analytical_position(Tend, 1.0, lambda)

            # Interface position error
            abs_pos_err = abs(xf_num - xf_exact_Tend)
            rel_pos_err = abs_pos_err / (abs(xf_exact_Tend) > 0 ? abs(xf_exact_Tend) : eps())

            # Extract Stefan residuals: max and last
            stefan_residual_max = NaN
            stefan_residual_last = NaN
            if !isempty(stefan_residuals)
                all_residuals = Float64[]
                for (_, resid_vec) in stefan_residuals
                    append!(all_residuals, resid_vec)
                end
                if !isempty(all_residuals)
                    stefan_residual_max = maximum(all_residuals)
                    stefan_residual_last = all_residuals[end]
                end
            end

            # Temperature errors (left of interface)
            u_num = solver.x[1:(nx + 1)]
            x_num = x_nodes_phys
            mask = x_num .<= xf_num

            if sum(mask) > 0
                u_anal_at_num = analytical_temperature.(x_num, Tend, 1.0, 1.0, lambda)
                u_num_below = u_num[mask]
                u_anal_below = u_anal_at_num[mask]
                L1_error = sum(abs.(u_num_below .- u_anal_below)) / length(u_num_below)
                L2_error = sqrt(sum((u_num_below .- u_anal_below).^2) / length(u_num_below))
            else
                L1_error = NaN
                L2_error = NaN
            end

            # Store results
            run_bucket = get!(results_by_tstart, Tstart, Vector{NamedTuple}())
            push!(run_bucket,
                (nx=nx, Tstart=Tstart, Tend=Tend, dx=dx, dt=dt,
                 L1_error=L1_error, L2_error=L2_error, rel_pos_error=rel_pos_err,
                 stefan_residual_max=stefan_residual_max, stefan_residual_last=stefan_residual_last))

            println("  -> L1=$L1_error, L2=$L2_error, rel_pos_err=$rel_pos_err, stefan_max=$stefan_residual_max OK")
        catch e
            println("  -> ERROR: $e")
            run_bucket = get!(results_by_tstart, Tstart, Vector{NamedTuple}())
            push!(run_bucket,
                (nx=nx, Tstart=Tstart, Tend=Tend, dx=dx, dt=dt,
                 L1_error=NaN, L2_error=NaN, rel_pos_error=NaN,
                 stefan_residual_max=NaN, stefan_residual_last=NaN))
        end
    end
end

# Save results to CSV, one file per Tstart with OOC columns
results_dir = joinpath(BENCH_ROOT, "results", "stefan")
mkpath(results_dir)

for (Tstart, runs) in sort(collect(results_by_tstart); by=first)
    df = DataFrame(runs)
    # Compute orders of convergence for L2 temperature and interface position
    ooc_L2 = compute_pairwise_orders(df.dx, df.L2_error)
    ooc_pos = compute_pairwise_orders(df.dx, df.rel_pos_error)
    df.ooc_L2 = ooc_L2
    df.ooc_pos = ooc_pos

    fname = "Stefan_1ph_convergence_Tstart$(Tstart).csv"
    csv_path = joinpath(results_dir, fname)
    CSV.write(csv_path, df)
    println("Saved: ", csv_path)
end

println("\nDone.")
