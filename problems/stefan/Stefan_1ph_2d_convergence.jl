using Penguin
using SpecialFunctions
using Roots
using Statistics
using CSV
using DataFrames

# Convergence study for a 2D one-phase Stefan test case with a height function.

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

# ============================================================================
# Analytical solution helpers
# ============================================================================

function find_lambda(Stefan_number)
    f = lambda -> sqrt(pi) * lambda * exp(lambda^2) * erf(lambda) - 1 / Stefan_number
    return find_zero(f, (1e-6, 10.0), Bisection())
end

function analytical_temperature(x, t, T0, k, lambda)
    t <= 0 && return T0
    return T0 - (T0 / erf(lambda)) * erf(x / (2 * sqrt(k * t)))
end

analytical_position(t, k, lambda) = 2 * lambda * sqrt(k * t)

# ============================================================================
# Convergence driver
# ============================================================================

function run_stefan_2d_convergence(nx_list::Vector{Int};
    lx::Float64=1.0,
    ly::Float64=1.0,
    x0::Float64=0.0,
    y0::Float64=0.0,
    Tstart::Float64=0.03,
    Tend::Float64=0.1,
    Stefan_number::Float64=1.0,
    T0::Float64=1.0,
    k::Float64=1.0,
    rho::Float64=1.0,
    L::Float64=1.0,
    max_iter::Int=3,
    tol::Float64=eps(),
    reltol::Float64=eps(),
    damping::Float64=1.0
)
    lambda = find_lambda(Stefan_number)

    h_vals = Float64[]
    dt_vals = Float64[]
    temp_L1 = Float64[]
    temp_L2 = Float64[]
    pos_rel = Float64[]
    stefan_max = Float64[]
    stefan_last = Float64[]

    for nx in nx_list
        ny = nx
        dx = lx / nx
        dt = 0.5 * dx^2

        mesh = Penguin.Mesh((nx, ny), (lx, ly), (x0, y0))
        x_offset = mesh.nodes[1][1][1] - x0
        STmesh = Penguin.SpaceTimeMesh(mesh, [Tstart, Tstart + dt], tag=mesh.tag)

        xf0_phys = analytical_position(Tstart, k, lambda)
        xf0 = xf0_phys + x_offset
        s_n(y) = xf0
        body = (x, y, t, _=0) -> (x - s_n(y))

        capacity = Capacity(body, STmesh; method="VOFI", integration_method=:vofijul)

        # Extract initial height information for the height-function interface
        Vn1 = capacity.A[3][1:end÷2, 1:end÷2]
        Vn = capacity.A[3][end÷2+1:end, end÷2+1:end]
        Vn = reshape(diag(Vn), (nx + 1, ny + 1))
        Vn1 = reshape(diag(Vn1), (nx + 1, ny + 1))

        Hn0 = collect(vec(sum(Vn, dims=1)))
        Hn1_0 = collect(vec(sum(Vn1, dims=1)))
        interface_position = x0 .+ Hn0 ./ (lx / nx)

        operator = DiffusionOps(capacity)
        bc = Dirichlet(0.0)
        bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(
            :top => Dirichlet(0.0),
            :bottom => Dirichlet(T0)
        ))
        stef_cond = InterfaceConditions(nothing, FluxJump(1.0, 1.0, rho * L))

        f = (x, y, z, t) -> 0.0
        K = (x, y, z) -> k
        fluid = Phase(capacity, operator, f, K)

        x_nodes = mesh.nodes[1][:, 1]
        x_nodes_phys = x_nodes .- x_offset

        u0_grid = zeros(nx + 1, ny + 1)
        for j in 1:(ny + 1)
            for i in 1:(nx + 1)
                u0_grid[i, j] = analytical_temperature(x_nodes_phys[i], Tstart, T0, k, lambda)
            end
        end

        u0o = vec(u0_grid)
        u0y = zeros((nx + 1) * (ny + 1))
        u0 = vcat(u0o, u0y)

        solver = MovingLiquidDiffusionUnsteadyMono(fluid, bc_b, bc, dt, u0, mesh, "BE")

        solver, residuals, xf_log, reconstruct, timestep_history = solve_MovingLiquidDiffusionUnsteadyMono2D!(
            solver, fluid, interface_position, Hn0, s_n, dt, Tstart, Tend, bc_b, bc, stef_cond, mesh, "BE";
            interpo="quad",
            Newton_params=(max_iter, tol, reltol, damping),
            adaptive_timestep=false,
            Δt_min=5e-4,
            method=Base.:\
        )

        xf_mean_log = [mean(xf_log[i]) for i in 1:length(xf_log)]
        Δt_eff = (Tend - Tstart) / length(xf_mean_log)
        times = Tstart .+ collect(1:length(xf_mean_log)) .* Δt_eff
        x_offset_0 = (xf_mean_log[1] - x_offset) - analytical_position(times[1], k, lambda)
        xf_log_phys = xf_mean_log .- x_offset .- x_offset_0
        xf_num = xf_log_phys[end]

        xf_exact = analytical_position(times[end], k, lambda)
        rel_pos_err = abs(xf_num - xf_exact) / (abs(xf_exact) > 0 ? abs(xf_exact) : eps())

        u_num = solver.x[1:(nx + 1) * (ny + 1)]
        y_center_idx = div(ny + 1, 2) + 1
        u_num_slice = u_num[(y_center_idx - 1) * (nx + 1) + 1:y_center_idx * (nx + 1)]
        mask = x_nodes_phys .<= xf_num

        if sum(mask) > 0
            u_anal = analytical_temperature.(x_nodes_phys, Tend, T0, k, lambda)
            u_num_below = u_num_slice[mask]
            u_anal_below = u_anal[mask]
            L1_error = sum(abs.(u_num_below .- u_anal_below)) / length(u_num_below)
            L2_error = sqrt(sum((u_num_below .- u_anal_below).^2) / length(u_num_below))
        else
            L1_error = NaN
            L2_error = NaN
        end

        stefan_residual_max = NaN
        stefan_residual_last = NaN
        if !isempty(residuals)
            all_residuals = Float64[]
            for (_, resid_vec) in residuals
                append!(all_residuals, resid_vec)
            end
            if !isempty(all_residuals)
                stefan_residual_max = maximum(all_residuals)
                stefan_residual_last = all_residuals[end]
            end
        end

        push!(h_vals, dx)
        push!(dt_vals, dt)
        push!(temp_L1, L1_error)
        push!(temp_L2, L2_error)
        push!(pos_rel, rel_pos_err)
        push!(stefan_max, stefan_residual_max)
        push!(stefan_last, stefan_residual_last)
    end

    ooc_L1 = compute_pairwise_orders(h_vals, temp_L1)
    ooc_L2 = compute_pairwise_orders(h_vals, temp_L2)
    ooc_pos = compute_pairwise_orders(h_vals, pos_rel)

    return (
        h_vals=h_vals,
        dt_vals=dt_vals,
        temp_L1=temp_L1,
        temp_L2=temp_L2,
        pos_rel=pos_rel,
        stefan_max=stefan_max,
        stefan_last=stefan_last,
        ooc_L1=ooc_L1,
        ooc_L2=ooc_L2,
        ooc_pos=ooc_pos
    )
end

function write_stefan_2d_csv(nx_list, data; csv_path=nothing)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "stefan") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "Stefan_1ph_2d_convergence.csv") : csv_path

    df = DataFrame(
        nx = nx_list,
        ny = nx_list,
        h = data.h_vals,
        dt = data.dt_vals,
        L1_error = data.temp_L1,
        L2_error = data.temp_L2,
        rel_pos_error = data.pos_rel,
        stefan_residual_max = data.stefan_max,
        stefan_residual_last = data.stefan_last,
        ooc_L1 = data.ooc_L1,
        ooc_L2 = data.ooc_L2,
        ooc_pos = data.ooc_pos
    )
    CSV.write(csv_out, df)
    return (csv_path=csv_out, table=df)
end

function main(; nx_list=nothing, csv_path=nothing)
    nx_vals = isnothing(nx_list) ? [4, 8, 16, 32] : nx_list
    data = run_stefan_2d_convergence(nx_vals)
    csv_info = write_stefan_2d_csv(nx_vals, data; csv_path=csv_path)
    return (data=data, csv_path=csv_info.csv_path, table=csv_info.table)
end

results = main()
println("Stefan 1ph 2d convergence written to:\n  ", results.csv_path)
