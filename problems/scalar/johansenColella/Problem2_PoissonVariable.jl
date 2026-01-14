using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
Johansen–Colella Problem 2: variable-coefficient Poisson equation on the same
star-shaped domain as Problem 1. The diffusion coefficient is β(r)=1−r² and the
equation ∇·(β ∇ϕ) = (7 r² − 15 r⁴) cos(3θ) retains the exact solution
ϕ(r, θ) = r⁴ cos(3θ). Interface Dirichlet data enforce the exact boundary values.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

center_default() = (0.5, 0.5)
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

function diffusivity_problem2(x::Float64, y::Float64, center=center_default())
    r2, _ = polar_features(x, y, center)
    return 1.0 - r2
end

function forcing_problem2(x::Float64, y::Float64, center=center_default())
    r2, θ = polar_features(x, y, center)
    return -(7.0 * r2 - 15.0 * r2^2) * cos(3θ)
end

function run_star_poisson_convergence(
    nx_list::Vector{Int},
    ny_list::Vector{Int},
    center::Tuple{Float64,Float64},
    u_exact::Function;
    source::Function,
    diffusivity::Function,
    lx::Float64 = 1.0,
    ly::Float64 = 1.0,
    relative::Bool = false
)
    h_vals = Float64[]
    err_vals_l1 = Float64[]
    err_full_l1 = Float64[]
    err_cut_l1 = Float64[]
    err_empty_l1 = Float64[]
    err_vals_linf = Float64[]
    err_full_linf = Float64[]
    err_cut_linf = Float64[]
    err_empty_linf = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()
    trunc_lp_all = Float64[]
    trunc_lp_full = Float64[]
    trunc_lp_cut = Float64[]
    trunc_linf_all = Float64[]
    trunc_linf_full = Float64[]
    trunc_linf_cut = Float64[]

    body = star_level_set

    for (nx, ny) in zip(nx_list, ny_list)
        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        capacity = Capacity(body, mesh; method="ImplicitIntegration")
        operator = DiffusionOps(capacity)

        bc_outer = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left   => bc_outer,
            :right  => bc_outer,
            :top    => bc_outer,
            :bottom => bc_outer
        ))
        phase = Phase(capacity, operator, source, diffusivity)
        solver = DiffusionSteadyMonoVariable(phase, bc_b, Dirichlet((x,y,_=0)->u_exact(x,y)))
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err_l1, full_err_l1, cut_err_l1, empty_err_l1 =
            check_convergence(u_exact, solver, capacity, 1, relative)
        _, _, global_err_linf, full_err_linf, cut_err_linf, empty_err_linf =
            check_convergence(u_exact, solver, capacity, Inf, relative)
        trunc = truncation_error_norms(solver, capacity, u_exact; norm=1)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals_l1, global_err_l1)
        push!(err_full_l1, full_err_l1)
        push!(err_cut_l1, cut_err_l1)
        push!(err_empty_l1, empty_err_l1)
        push!(err_vals_linf, global_err_linf)
        push!(err_full_linf, full_err_linf)
        push!(err_cut_linf, cut_err_linf)
        push!(err_empty_linf, empty_err_linf)
        push!(inside_cells, count_inside_cells(capacity))
        Δx = lx / nx
        Δy = ly / ny
        coverage_x = ceil(Int, 2 * star_max_radius() / Δx)
        coverage_y = ceil(Int, 2 * star_max_radius() / Δy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
        push!(trunc_lp_all, trunc.lp_all)
        push!(trunc_lp_full, trunc.lp_full)
        push!(trunc_lp_cut, trunc.lp_cut)
        push!(trunc_linf_all, trunc.linf_all)
        push!(trunc_linf_full, trunc.linf_full)
        push!(trunc_linf_cut, trunc.linf_cut)
    end

    return (
        l1 = (
            h_vals = h_vals,
            err_vals = err_vals_l1,
            err_full_vals = err_full_l1,
            err_cut_vals = err_cut_l1,
            err_empty_vals = err_empty_l1,
            inside_cells = inside_cells,
            inside_cells_by_dim = inside_cells_by_dim,
            trunc_lp_all = trunc_lp_all,
            trunc_lp_full = trunc_lp_full,
            trunc_lp_cut = trunc_lp_cut,
            trunc_linf_all = trunc_linf_all,
            trunc_linf_full = trunc_linf_full,
            trunc_linf_cut = trunc_linf_cut,
            orders = compute_orders(h_vals, err_vals_l1, err_full_l1, err_cut_l1),
            norm = 1
        ),
        linf = (
            h_vals = h_vals,
            err_vals = err_vals_linf,
            err_full_vals = err_full_linf,
            err_cut_vals = err_cut_linf,
            err_empty_vals = err_empty_linf,
            inside_cells = inside_cells,
            inside_cells_by_dim = inside_cells_by_dim,
            trunc_lp_all = trunc_lp_all,
            trunc_lp_full = trunc_lp_full,
            trunc_lp_cut = trunc_lp_cut,
            trunc_linf_all = trunc_linf_all,
            trunc_linf_full = trunc_linf_full,
            trunc_linf_cut = trunc_linf_cut,
            orders = compute_orders(h_vals, err_vals_linf, err_full_linf, err_cut_linf),
            norm = Inf
        )
    )
end

function write_convergence_csv(method_name, data_l1, data_linf; csv_path=nothing)
    df_l1 = make_convergence_dataframe(method_name, data_l1)
    df_linf = make_convergence_dataframe(method_name, data_linf)
    df = vcat(df_l1, df_linf)
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
        source = (x,y,_)->forcing_problem2(x,y,center),
        diffusivity = (x,y,_)->diffusivity_problem2(x,y,center),
        lx = 1.0,
        ly = 1.0
    )

    csv_info = write_convergence_csv("JohansenColella_P2", data.l1, data.linf; csv_path=csv_path)
    return (data_l1 = data.l1, data_linf = data.linf, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

# Plot

@testset "Johansen-Colella Problem 2" begin
    orders = results.data_l1.orders
    @test !isnan(orders.all)
    @test orders.all > 1.0
    @test length(results.data_l1.h_vals) == length(results.data_l1.err_vals)
    @test results.data_l1.h_vals[1] > results.data_l1.h_vals[end]
    @test minimum(results.data_l1.err_vals) < maximum(results.data_l1.err_vals)
    @test isfile(results.csv_path)
    println(orders)
end

function read_johansen_colella_table3(table_path::AbstractString)
    n_vals = Int[]
    e_max = Float64[]
    e_avg = Float64[]
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
        push!(e_max, parse(Float64, fields[3]))
        push!(e_avg, parse(Float64, fields[4]))
    end
    return (n_vals = n_vals, e_max = e_max, e_avg = e_avg)
end

using CairoMakie
function plot_solution_convergence_with_johansen(;
    data_l1 = results.data_l1,
    data_linf = results.data_linf,
    table_path = joinpath(BENCH_ROOT, "results", "scalar", "johansencolella", "neumann.table3"),
    output_path = "johansen_colella_problem2_convergence_comparison.png"
)
    jc = read_johansen_colella_table3(table_path)
    h_jc = 1.0 ./ jc.n_vals

    h_vals = data_l1.h_vals
    fig = Figure(resolution=(640, 480))
    ax = Axis(
        fig[1, 1],
        title = "Log-Log Solution Error Convergence (Problem 2)",
        xlabel = "h",
        ylabel = "Error norm",
        xscale = log10,
        yscale = log10
    )

    lines!(ax, h_vals, data_l1.err_vals, label = "All cells L1 (ours)")
    scatter!(ax, h_vals, data_l1.err_vals, marker = :circle)
    lines!(ax, h_vals, data_linf.err_vals, label = "All cells Linf (ours)")
    scatter!(ax, h_vals, data_linf.err_vals, marker = :square)

    lines!(ax, h_jc, jc.e_avg, linestyle = :dash, label = "JC e_avg (L1)")
    scatter!(ax, h_jc, jc.e_avg, marker = :utriangle)
    lines!(ax, h_jc, jc.e_max, linestyle = :dash, label = "JC e_max (Linf)")
    scatter!(ax, h_jc, jc.e_max, marker = :diamond)

    axislegend(ax, position = :rb)
    display(fig)
    save(output_path, fig)
end

plot_solution_convergence_with_johansen()

function plot_solutionn(solver, mesh, capacity; center=center_default())
    u_ana, u_num, _, _, _, _ = check_convergence((x,y)->exact_phi(x, y, center), solver, capacity)

    cell_types = capacity.cell_types
    u_num_promasked = copy(u_num)
    u_ana_promasked = copy(u_ana)
    u_num_promasked[cell_types .== 0] .= NaN
    u_ana_promasked[cell_types .== 0] .= NaN

    xg = mesh.nodes[1]
    yg = mesh.nodes[2]
    nx_plot = length(xg)
    ny_plot = length(yg)

    sol_grid = reshape(u_num_promasked, (nx_plot, ny_plot))'
    ana_grid = reshape(u_ana_promasked, (nx_plot, ny_plot))'
    err_grid = sol_grid .- ana_grid

    fig_sol = Figure(resolution=(500,400))
    ax1 = Axis(fig_sol[1,1], title="Solution (numerical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm1 = heatmap!(ax1, xg, yg, sol_grid, colormap = :viridis)
    Colorbar(fig_sol[1,2], hm1, label = "Value")
    display(fig_sol)
    save("johansen_colella_problem2_solution.png", fig_sol)

    fig_err = Figure(resolution=(500,400))
    ax3 = Axis(fig_err[1,1], title="Error (numerical - analytical)", xlabel="x", ylabel="y", aspect=DataAspect())
    hm3 = heatmap!(ax3, xg, yg, err_grid, colormap = :viridis)
    Colorbar(fig_err[1,2], hm3, label = "Error")
    display(fig_err)

    fig_logerr = Figure(resolution=(500,400))
    ax4 = Axis(fig_logerr[1,1], title="Log10 Error Magnitude", xlabel="x", ylabel="y", aspect=DataAspect())
    log_err_grid = log10.(abs.(err_grid))
    hm4 = heatmap!(ax4, xg, yg, log_err_grid, colormap = :viridis)
    Colorbar(fig_logerr[1,2], hm4, label = "Log10 |Error|")
    display(fig_logerr)
    save("johansen_colella_problem2_log_error.png", fig_logerr)
end

function plot_solution_and_error(; nx=640, ny=640, center=center_default())
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
        (x,y,_)->forcing_problem2(x, y, center),
        (x,y,_)->diffusivity_problem2(x, y, center)
    )
    solver = DiffusionSteadyMonoVariable(phase, bc_b, Dirichlet((x,y,_=0)->exact_phi(x, y, center)))
    solve_DiffusionSteadyMono!(solver; method=Base.:\)

    plot_solutionn(solver, mesh, capacity; center=center)
end

plot_solution_and_error()
