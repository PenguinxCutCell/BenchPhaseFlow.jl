using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
Johansen-Colella Problem 1 (Neumann cut): constant-coefficient Poisson equation
outside a star-shaped boundary defined by
    Omega = {(r, theta) : r >= 0.25 + 0.05 cos(6 theta)}.
We solve Delta phi = 7 r^2 cos(3 theta) with Dirichlet data taken from the exact
solution phi(r, theta) = r^4 cos(3 theta) on the outer box, and Neumann data on
 the cut boundary from the analytical gradient.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

center_default() = (0.49, 0.5)
base_radius() = 0.25
perturb_radius() = 0.05
star_max_radius() = base_radius() + perturb_radius()

function star_level_set(x, y, _=0)
    dx = x - center_default()[1]
    dy = y - center_default()[2]
    theta = atan(dy, dx)
    r = sqrt(dx^2 + dy^2)
    return (base_radius() + perturb_radius() * cos(6 * theta)) - r
end

function polar_features(x::Float64, y::Float64, center::Tuple{Float64,Float64})
    dx = x - center[1]
    dy = y - center[2]
    r = sqrt(dx^2 + dy^2)
    theta = atan(dy, dx)
    return r, theta
end

function exact_phi(x::Float64, y::Float64, center=center_default())
    r, theta = polar_features(x, y, center)
    return (r^4) * cos(3 * theta)
end

function forcing_problem1(x::Float64, y::Float64, center=center_default())
    r, theta = polar_features(x, y, center)
    return -7.0 * (r^2) * cos(3 * theta)
end

function star_normal_from_level_set(x::Float64, y::Float64, center::Tuple{Float64,Float64})
    r, theta = polar_features(x, y, center)
    dphidr = -1.0
    dphidtheta_over_r = -6.0 * perturb_radius() * sin(6 * theta) / r
    nx = dphidr * cos(theta) - dphidtheta_over_r * sin(theta)
    ny = dphidr * sin(theta) + dphidtheta_over_r * cos(theta)
    norm = sqrt(nx^2 + ny^2)
    return nx / norm, ny / norm
end

function normal_from_args(x::Float64, y::Float64, center::Tuple{Float64,Float64}, args...)
    if length(args) == 1
        arg = args[1]
        if arg isa Tuple && length(arg) >= 2
            return arg[1], arg[2]
        elseif arg isa AbstractVector && length(arg) >= 2
            return arg[1], arg[2]
        end
    elseif length(args) >= 2
        return args[end - 1], args[end]
    end
    return star_normal_from_level_set(x, y, center)
end

function exact_normal_flux(x::Float64, y::Float64, center::Tuple{Float64,Float64}, args...)
    r, theta = polar_features(x, y, center)
    dsdr = 4.0 * r^3 * cos(3 * theta)
    dsdtheta = -3.0 * r^3 * sin(3 * theta)

    nx, ny = normal_from_args(x, y, center, args...)
    norm = sqrt(nx^2 + ny^2)
    nx /= norm
    ny /= norm

    return nx * (dsdr * cos(theta) - dsdtheta * sin(theta)) +
           ny * (dsdr * sin(theta) + dsdtheta * cos(theta))
end

function run_star_poisson_convergence(
    nx_list::Vector{Int},
    ny_list::Vector{Int},
    center::Tuple{Float64,Float64},
    u_exact::Function;
    source::Function,
    diffusivity::Function = (x,y,_)->1.0,
    lx::Float64 = 1.0,
    ly::Float64 = 1.0,
    norm::Int = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full = Float64[]
    err_cut = Float64[]
    err_empty = Float64[]
    trunc_l2_all = Float64[]
    trunc_l2_full = Float64[]
    trunc_l2_cut = Float64[]
    trunc_linf_all = Float64[]
    trunc_linf_full = Float64[]
    trunc_linf_cut = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    for (nx, ny) in zip(nx_list, ny_list)
        dx = lx / nx
        dy = ly / ny
        x_shift = dx
        y_shift = dy
        Lx_eff = lx - dx
        Ly_eff = ly - dy

        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        body = (x, y, _=0) -> star_level_set(x - x_shift, y - y_shift)
        capacity = Capacity(body, mesh; method="VOFI", integration_method=:vofijul)
        operator = DiffusionOps(capacity)
        u_exact_shifted = (x, y) -> u_exact((x - x_shift)/Lx_eff, (y - y_shift)/Ly_eff)
        source_shifted = (x, y, _=0) -> source((x - x_shift)/Lx_eff, (y - y_shift)/Ly_eff)
        diffusivity_shifted = (x, y, _=0) -> diffusivity((x - x_shift)/Lx_eff, (y - y_shift)/Ly_eff)

        bc_outer = Dirichlet((x,y,_=0)->u_exact_shifted(x,y))
        bc_b = BorderConditions(Dict(
            :left   => bc_outer,
            :right  => bc_outer,
            :top    => bc_outer,
            :bottom => bc_outer
        ))
        phase = Phase(capacity, operator, source_shifted, diffusivity_shifted)
        bc_interface = Neumann((x,y,args...)->exact_normal_flux(x - x_shift, y - y_shift, center, args...))
        solver = DiffusionSteadyMono(phase, bc_b, bc_interface)
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact_shifted, solver, capacity, norm, relative)
        trunc = truncation_error_norms(solver, capacity, u_exact_shifted)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, global_err)
        push!(err_full, full_err)
        push!(err_cut, cut_err)
        push!(err_empty, empty_err)
        push!(trunc_l2_all, trunc.l2_all)
        push!(trunc_l2_full, trunc.l2_full)
        push!(trunc_l2_cut, trunc.l2_cut)
        push!(trunc_linf_all, trunc.linf_all)
        push!(trunc_linf_full, trunc.linf_full)
        push!(trunc_linf_cut, trunc.linf_cut)
        push!(inside_cells, count_inside_cells(capacity))
        coverage_x = ceil(Int, 2 * star_max_radius() / dx)
        coverage_y = ceil(Int, 2 * star_max_radius() / dy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full,
        err_cut_vals = err_cut,
        err_empty_vals = err_empty,
        trunc_l2_all = trunc_l2_all,
        trunc_l2_full = trunc_l2_full,
        trunc_l2_cut = trunc_l2_cut,
        trunc_linf_all = trunc_linf_all,
        trunc_linf_full = trunc_linf_full,
        trunc_linf_cut = trunc_linf_cut,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        orders = compute_orders(h_vals, err_vals, err_full, err_cut),
        trunc_orders_l2 = compute_orders(h_vals, trunc_l2_all, trunc_l2_full, trunc_l2_cut),
        trunc_orders_linf = compute_orders(h_vals, trunc_linf_all, trunc_linf_full, trunc_linf_cut),
        norm = norm
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_convergence_dataframe(method_name, data)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "scalar") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "$(method_name)_Convergence.csv") : csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing, ny_list=nothing)
    center = center_default()
    nx_vals = isnothing(nx_list) ? [16, 32, 64, 128, 256, 512] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list
    u_exact = (x,y) -> exact_phi(x, y, center)

    data = run_star_poisson_convergence(
        nx_vals, ny_vals, center, u_exact;
        source = (x,y,_=0)->forcing_problem1(x, y, center),
        diffusivity = (x,y,_=0)->1.0,
        lx = 1.0,
        ly = 1.0,
        norm = 1
    )

    csv_info = write_convergence_csv("JohansenColella_P1_Neumann", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Johansen-Colella Problem 1 Neumann" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    #@test orders.all > 1.0
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
    println(orders)
end


# Using cairomakie, plot the solution and the error on the finest mesh
using CairoMakie
function plot_solutionn(solver, mesh, capacity; center=center_default(), shift=(0.0, 0.0))

    # Use check_convergence to obtain analytical and numerical cell-centered fields
    u_exact = (x, y) -> exact_phi(x - shift[1], y - shift[2], center)
    u_ana, u_num, _, _, _, _ = check_convergence(u_exact, solver, capacity)

    # Mask inactive cells (empty cells) using capacity.cell_types
    cell_types = capacity.cell_types
    u_num_promasked = copy(u_num)
    u_ana_promasked = copy(u_ana)
    u_num_promasked[cell_types .== 0] .= NaN
    u_ana_promasked[cell_types .== 0] .= NaN

    # Prepare plotting grid and reshape arrays to (ny, nx) as used by heatmap
    xg = mesh.nodes[1]
    yg = mesh.nodes[2]
    nx_plot = length(xg)
    ny_plot = length(yg)


    sol_grid = reshape(u_num_promasked, (nx_plot, ny_plot))'
    ana_grid = reshape(u_ana_promasked, (nx_plot, ny_plot))'
    err_grid = sol_grid .- ana_grid

    # Solution figure
    fig_sol = Figure(resolution=(500,400))
    ax1 = Axis(fig_sol[1,1], title="Solution (numerical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm1 = heatmap!(ax1, xg, yg, sol_grid, colormap = :viridis)
    Colorbar(fig_sol[1,2], hm1, label = "Value")
    display(fig_sol)
    save("johansen_colella_problem1_neumann_solution.png", fig_sol)

    # Analytical figure (separate)
    fig_ana = Figure(resolution=(500,400))
    ax2 = Axis(fig_ana[1,1], title="Solution (analytical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm2 = heatmap!(ax2, xg, yg, ana_grid, colormap = :viridis)
    Colorbar(fig_ana[1,2], hm2, label = "Value")
    display(fig_ana)
    save("johansen_colella_problem1_neumann_analytical.png", fig_ana)

    # Error figure (separate)
    fig_err = Figure(resolution=(500,400))
    ax3 = Axis(fig_err[1,1], title="Error (numerical - analytical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm3 = heatmap!(ax3, xg, yg, err_grid, colormap = :viridis)
    Colorbar(fig_err[1,2], hm3, label = "Error")
    display(fig_err)

    # Log error figure (separate)
    fig_logerr = Figure(resolution=(500,400))
    ax4 = Axis(fig_logerr[1,1], title="Log10 Error Magnitude", xlabel="x", ylabel="y", aspect=DataAspect())
    log_err_grid = log10.(abs.(err_grid))
    hm4 = heatmap!(ax4, xg, yg, log_err_grid, colormap = :viridis)
    Colorbar(fig_logerr[1,2], hm4, label = "Log10 |Error|")
    display(fig_logerr)
    save("johansen_colella_problem1_neumann_log_error.png", fig_logerr)
end

function plot_solution_and_error(; nx=512, ny=512, center=center_default())
    lx = 1.0
    ly = 1.0
    dx = lx / nx
    dy = ly / ny
    shift = (dx, dy)

    mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
    body = (x, y, _=0) -> star_level_set(x - shift[1], y - shift[2])
    capacity = Capacity(body, mesh; method="ImplicitIntegration")
    operator = DiffusionOps(capacity)
    bc_outer = Dirichlet((x,y,_=0)->exact_phi(x - shift[1], y - shift[2], center))
    bc_b = BorderConditions(Dict(
        :left   => bc_outer,
        :right  => bc_outer,
        :top    => bc_outer,
        :bottom => bc_outer
    ))
    phase = Phase(capacity, operator,
        (x,y,_)->forcing_problem1(x - shift[1], y - shift[2], center),
        (x,y,_)->1.0
    )
    bc_interface = Neumann((x,y,args...)->exact_normal_flux(x - shift[1], y - shift[2], center, args...))
    solver = DiffusionSteadyMono(phase, bc_b, bc_interface)
    solve_DiffusionSteadyMono!(solver; method=Base.:\)

    plot_solutionn(solver, mesh, capacity; center=center, shift=shift)
end

plot_solution_and_error()
