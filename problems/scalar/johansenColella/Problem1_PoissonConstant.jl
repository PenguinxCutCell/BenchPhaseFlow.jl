using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using DataFrames
using Test

"""
Johansen–Colella Problem 1: constant-coefficient Poisson equation inside a
star-shaped domain defined by
    Ω = {(r, θ) : r ≤ 0.30 + 0.15 cos(6θ)}.
We solve Δϕ = 7 r² cos(3θ) with Dirichlet data taken from the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Convergence is measured on uniform Cartesian meshes using
Penguin's implicit-integration capacity.

Truncation errors are computed by inserting the exact solution sampled at the cut-cell centroids into the discrete operator and comparing it with the continuous operator evaluated at the same centroid locations. 
This differs slightly from the convention used by Johansen and Colella, who evaluate the exact solution at Cartesian cell centers while evaluating the right-hand side at cut-cell centroids. 
The error constants might differ.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

center_default() = (0.49, 0.5)
base_radius() = 0.30
perturb_radius() = 0.15
star_max_radius() = base_radius() + perturb_radius()

function star_level_set(x, y, _=0)
    dx = x - center_default()[1]
    dy = y - center_default()[2]
    θ = atan(dy, dx)
    r = sqrt(dx^2 + dy^2)
    return r - (base_radius() + perturb_radius() * cos(6θ))
end

function polar_features(x::Float64, y::Float64, center::Tuple{Float64,Float64})
    dx = x - center[1]
    dy = y - center[2]
    r2 = dx^2 + dy^2
    θ = atan(dy, dx)
    return r2, θ
end

function exact_phi(x::Float64, y::Float64, center=center_default())
    r2, θ = polar_features(x, y, center)
    return (r2^2) * cos(3θ)
end

function forcing_problem1(x::Float64, y::Float64, center=center_default())
    r2, θ = polar_features(x, y, center)
    return -7.0 * r2 * cos(3θ)
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
    trunc_lp_all = Float64[]
    trunc_lp_full = Float64[]
    trunc_lp_cut = Float64[]
    trunc_linf_all = Float64[]
    trunc_linf_full = Float64[]
    trunc_linf_cut = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    body = star_level_set

    for (nx, ny) in zip(nx_list, ny_list)
        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        capacity = Capacity(body, mesh; method="VOFI",integration_method=:vofijul)
        operator = DiffusionOps(capacity)

        bc_outer = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left   => bc_outer,
            :right  => bc_outer,
            :top    => bc_outer,
            :bottom => bc_outer
        ))
        phase = Phase(capacity, operator, source, diffusivity)
        solver = DiffusionSteadyMono(phase, bc_b, Dirichlet((x,y,_=0)->u_exact(x,y)))
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)
        trunc = truncation_error_norms(solver, capacity, u_exact)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, global_err)
        push!(err_full, full_err)
        push!(err_cut, cut_err)
        push!(err_empty, empty_err)
        push!(trunc_lp_all, trunc.lp_all)
        push!(trunc_lp_full, trunc.lp_full)
        push!(trunc_lp_cut, trunc.lp_cut)
        push!(trunc_linf_all, trunc.linf_all)
        push!(trunc_linf_full, trunc.linf_full)
        push!(trunc_linf_cut, trunc.linf_cut)
        push!(inside_cells, count_inside_cells(capacity))
        Δx = lx / nx
        Δy = ly / ny
        coverage_x = ceil(Int, 2 * star_max_radius() / Δx)
        coverage_y = ceil(Int, 2 * star_max_radius() / Δy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full,
        err_cut_vals = err_cut,
        err_empty_vals = err_empty,
        trunc_lp_all = trunc_lp_all,
        trunc_lp_full = trunc_lp_full,
        trunc_lp_cut = trunc_lp_cut,
        trunc_linf_all = trunc_linf_all,
        trunc_linf_full = trunc_linf_full,
        trunc_linf_cut = trunc_linf_cut,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        orders = compute_orders(h_vals, err_vals, err_full, err_cut),
        trunc_orders_lp = compute_orders(h_vals, trunc_lp_all, trunc_lp_full, trunc_lp_cut),
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
    nx_vals = isnothing(nx_list) ? [40, 80, 160, 320, 640] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list
    u_exact = (x,y) -> exact_phi(x,y,center)

    data = run_star_poisson_convergence(
        nx_vals, ny_vals, center, u_exact;
        source = (x,y,_)->forcing_problem1(x,y,center),
        diffusivity = (x,y,_)->1.0,
        lx = 1.0,
        ly = 1.0,
        norm = 1
    )

    csv_info = write_convergence_csv("JohansenColella_P1", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

"""
@testset "Johansen-Colella Problem 1" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test orders.all > 1.0
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
    println(orders)
end
"""

function read_johansen_colella_table(table_path::AbstractString)
    n_vals = Int[]
    e_max = Float64[]
    e_p_avg = Float64[]
    e_f_avg = Float64[]
    for line in eachline(table_path)
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, "#")
            continue
        end
        fields = split(stripped)
        if length(fields) < 4
            continue
        end
        push!(n_vals, parse(Int, fields[1]))
        push!(e_max, parse(Float64, fields[2]))
        push!(e_p_avg, parse(Float64, fields[3]))
        push!(e_f_avg, parse(Float64, fields[4]))
    end
    return (n_vals = n_vals, e_max = e_max, e_p_avg = e_p_avg, e_f_avg = e_f_avg)
end

function read_convergence_csv(csv_path::AbstractString)
    df = CSV.read(csv_path, DataFrame)
    h_vals = collect(df.h)
    err_vals = collect(df.all_err)
    err_full_vals = collect(df.full_err)
    err_cut_vals = collect(df.cut_err)
    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full_vals,
        err_cut_vals = err_cut_vals
    )
end

using CairoMakie
function plot_convergence_comparison(;
    csv_path = joinpath(BENCH_ROOT, "results", "scalar", "JohansenColella_P1_Convergence.csv"),
    table_path = joinpath(BENCH_ROOT, "results", "scalar", "johansencolella", "neumann.table2"),
    output_path = "johansen_colella_problem1_convergence_comparison.png"
)
    results_data = read_convergence_csv(csv_path)
    jc = read_johansen_colella_table(table_path)
    h_jc = 1.0 ./ jc.n_vals
    h_vals = results_data.h_vals

    fig = Figure(resolution=(640, 480))
    ax = Axis(
        fig[1, 1],
        title = "Log-Log Convergence: Problem 1 vs Johansen-Colella",
        xlabel = "h",
        ylabel = "Error norm",
        xscale = log10,
        yscale = log10
    )

    lines!(ax, h_vals, results_data.err_vals,  label = "All cells (ours)")
    scatter!(ax, h_vals, results_data.err_vals, marker = :circle)
    lines!(ax, h_vals, results_data.err_full_vals,  label = "Full cells (ours)")
    scatter!(ax, h_vals, results_data.err_full_vals, marker = :square)
    lines!(ax, h_vals, results_data.err_cut_vals, label = "Cut cells (ours)")
    scatter!(ax, h_vals, results_data.err_cut_vals, marker = :diamond)

    lines!(ax, h_jc, jc.e_max, linestyle = :dash, label = "JC e_max")
    scatter!(ax, h_jc, jc.e_max, marker = :utriangle)
    lines!(ax, h_jc, jc.e_p_avg, linestyle = :dash, label = "JC e_p_avg")
    scatter!(ax, h_jc, jc.e_p_avg, marker = :hexagon)
    lines!(ax, h_jc, jc.e_f_avg, linestyle = :dash, label = "JC e_f_avg")
    scatter!(ax, h_jc, jc.e_f_avg, marker = :star5)

    # Reference slopes
    x_ref = [minimum(h_vals)/2, maximum(h_vals)*2]
    y_ref_2 = 0.1 .* x_ref .^ 2
    y_ref_3 = 0.1 .* x_ref .^ 3
    lines!(ax, x_ref, y_ref_2, color = :black, linestyle = :dot, label = "O(h²)")
    lines!(ax, x_ref, y_ref_3, color = :black, linestyle = :dashdot, label = "O(h³)")

    axislegend(ax, position = :rb)
    display(fig)
    save(output_path, fig)
end

plot_convergence_comparison()

function read_johansen_colella_trunc_table(table_path::AbstractString)
    n_vals = Int[]
    tau_p_max = Float64[]
    tau_p_avg = Float64[]
    tau_f_max = Float64[]
    for line in eachline(table_path)
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, "#")
            continue
        end
        fields = split(stripped)
        if length(fields) < 4
            continue
        end
        push!(n_vals, parse(Int, fields[1]))
        push!(tau_p_max, parse(Float64, fields[2]))
        push!(tau_p_avg, parse(Float64, fields[3]))
        push!(tau_f_max, parse(Float64, fields[4]))
    end
    return (n_vals = n_vals, tau_p_max = tau_p_max, tau_p_avg = tau_p_avg, tau_f_max = tau_f_max)
end

function read_truncation_convergence_csv(csv_path::AbstractString)
    df = CSV.read(csv_path, DataFrame)
    return (
        h_vals = collect(df.h),
        trunc_all = collect(df.trunc_all),
        trunc_full = collect(df.trunc_full),
        trunc_cut = collect(df.trunc_cut),
        trunc_max_all = collect(df.trunc_max_all),
        trunc_max_full = collect(df.trunc_max_full),
        trunc_max_cut = collect(df.trunc_max_cut)
    )
end

function plot_truncation_convergence_comparison(;
    csv_path = joinpath(BENCH_ROOT, "results", "scalar", "JohansenColella_P1_Convergence.csv"),
    table_path = joinpath(BENCH_ROOT, "results", "scalar", "johansencolella", "neumann.table1"),
    output_path = "johansen_colella_problem1_truncation_convergence.png"
)
    trunc = read_truncation_convergence_csv(csv_path)
    jc = read_johansen_colella_trunc_table(table_path)
    h_jc = 1.0 ./ jc.n_vals

    fig = Figure(resolution=(640, 480))
    ax = Axis(
        fig[1, 1],
        title = "Log-Log Truncation Error Convergence",
        xlabel = "h",
        ylabel = "Truncation error norm",
        xscale = log10,
        yscale = log10
    )

    lines!(ax, trunc.h_vals, trunc.trunc_all, label = "Lp all (ours)")
    scatter!(ax, trunc.h_vals, trunc.trunc_all, marker = :circle)
    lines!(ax, trunc.h_vals, trunc.trunc_full, label = "Lp full (ours)")
    scatter!(ax, trunc.h_vals, trunc.trunc_full, marker = :square)
    lines!(ax, trunc.h_vals, trunc.trunc_cut, label = "Lp cut (ours)")
    scatter!(ax, trunc.h_vals, trunc.trunc_cut, marker = :diamond)

    lines!(ax, trunc.h_vals, trunc.trunc_max_all, label = "L∞ all (ours)")
    scatter!(ax, trunc.h_vals, trunc.trunc_max_all, marker = :utriangle)
    lines!(ax, trunc.h_vals, trunc.trunc_max_full, label = "L∞ full (ours)")
    scatter!(ax, trunc.h_vals, trunc.trunc_max_full, marker = :dtriangle)
    lines!(ax, trunc.h_vals, trunc.trunc_max_cut, label = "L∞ cut (ours)")
    scatter!(ax, trunc.h_vals, trunc.trunc_max_cut, marker = :pentagon)

    lines!(ax, h_jc, jc.tau_p_max, linestyle = :dash, label = "JC tau_p_max")
    scatter!(ax, h_jc, jc.tau_p_max, marker = :hexagon)
    lines!(ax, h_jc, jc.tau_p_avg, linestyle = :dash, label = "JC tau_p_avg")
    scatter!(ax, h_jc, jc.tau_p_avg, marker = :star5)
    lines!(ax, h_jc, jc.tau_f_max, linestyle = :dash, label = "JC tau_f_max")
    scatter!(ax, h_jc, jc.tau_f_max, marker = :cross)

    axislegend(ax, position = :rb)
    display(fig)
    save(output_path, fig)
end

plot_truncation_convergence_comparison()

function read_basilisk_neumann_ref(table_path::AbstractString)
    n_vals = Int[]
    n_avg = Float64[]
    n_max = Float64[]
    np_avg = Float64[]
    np_max = Float64[]
    nf_avg = Float64[]
    nf_max = Float64[]
    for line in eachline(table_path)
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, "#") || startswith(stripped, "maxres")
            continue
        end
        fields = split(stripped)
        if length(fields) < 7
            continue
        end
        push!(n_vals, parse(Int, fields[1]))
        push!(n_avg, parse(Float64, fields[2]))
        push!(n_max, parse(Float64, fields[3]))
        push!(np_avg, parse(Float64, fields[4]))
        push!(np_max, parse(Float64, fields[5]))
        push!(nf_avg, parse(Float64, fields[6]))
        push!(nf_max, parse(Float64, fields[7]))
    end
    return (
        n_vals = n_vals,
        n_avg = n_avg,
        n_max = n_max,
        np_avg = np_avg,
        np_max = np_max,
        nf_avg = nf_avg,
        nf_max = nf_max
    )
end

function plot_solution_convergence_with_basilisk(;
    csv_path = joinpath(BENCH_ROOT, "results", "scalar", "JohansenColella_P1_Convergence.csv"),
    table_path = joinpath(BENCH_ROOT, "results", "scalar", "johansencolella", "neumann.ref"),
    output_path = "johansen_colella_problem1_basilisk_convergence.png"
)
    ours = read_convergence_csv(csv_path)
    bas = read_basilisk_neumann_ref(table_path)
    h_bas = 1.0 ./ bas.n_vals

    fig = Figure(resolution=(640, 480))
    ax = Axis(
        fig[1, 1],
        title = "Log-Log Solution Error Convergence (Basilisk comparison)",
        xlabel = "h",
        ylabel = "Error norm",
        xscale = log10,
        yscale = log10
    )

    lines!(ax, ours.h_vals, ours.err_vals, label = "All cells L2 (ours)")
    scatter!(ax, ours.h_vals, ours.err_vals, marker = :circle)
    lines!(ax, ours.h_vals, ours.err_full_vals, label = "Full cells L2 (ours)")
    scatter!(ax, ours.h_vals, ours.err_full_vals, marker = :square)
    lines!(ax, ours.h_vals, ours.err_cut_vals, label = "Cut cells L2 (ours)")
    scatter!(ax, ours.h_vals, ours.err_cut_vals, marker = :diamond)

    lines!(ax, h_bas, bas.n_avg, linestyle = :dash, label = "Basilisk total avg")
    scatter!(ax, h_bas, bas.n_avg, marker = :utriangle)
    lines!(ax, h_bas, bas.np_avg, linestyle = :dash, label = "Basilisk partial avg")
    scatter!(ax, h_bas, bas.np_avg, marker = :hexagon)
    lines!(ax, h_bas, bas.nf_avg, linestyle = :dash, label = "Basilisk full avg")
    scatter!(ax, h_bas, bas.nf_avg, marker = :star5)

    axislegend(ax, position = :rb)
    display(fig)
    save(output_path, fig)
end

plot_solution_convergence_with_basilisk()

# Using cairomakie, plot the solution and the error on the finest mesh
function plot_solutionn(solver, mesh, body::Function, capacity; state_i=1)
   
    # Use check_convergence to obtain analytical and numerical cell-centered fields
    u_ana, u_num, _, _, _, _ = check_convergence((x,y)->exact_phi(x,y), solver, capacity)

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
    save("johansen_colella_problem1_solution.png", fig_sol)

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
    save("johansen_colella_problem1_log_error.png", fig_logerr)
end
function plot_solution_and_error(; nx=512, ny=512, center=center_default())
    lx = 1.0
    ly = 1.0
    mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
    capacity = Capacity(star_level_set, mesh; method="ImplicitIntegration")
    operator = DiffusionOps(capacity) 
    bc_outer = Dirichlet(0.0)
    bc_b = BorderConditions(Dict(
        :left   => bc_outer,
        :right  => bc_outer,
        :top    => bc_outer,
        :bottom => bc_outer
    ))
    phase = Phase(capacity, operator,
        (x,y,_)->forcing_problem1(x,y,center),
        (x,y,_)->1.0
    )
    solver = DiffusionSteadyMono(phase, bc_b, Dirichlet((x,y,_=0)->exact_phi(x,y,center)))
    solve_DiffusionSteadyMono!(solver; method=Base.:\)
    
    plot_solutionn(solver, mesh, star_level_set, capacity)
end

plot_solution_and_error()
