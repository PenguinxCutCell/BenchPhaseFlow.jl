using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using DataFrames
using Test

"""
Johansen-Colella Problem 2: truncation error comparison for mean types.
Compares L1 and Linf truncation errors for default and arithmetic means.
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
    theta = atan(dy, dx)
    r = sqrt(dx^2 + dy^2)
    return r - (base_radius() + perturb_radius() * cos(6 * theta))
end

function polar_features(x::Float64, y::Float64, center::Tuple{Float64,Float64})
    dx = x - center[1]
    dy = y - center[2]
    r2 = dx^2 + dy^2
    theta = atan(dy, dx)
    return r2, theta
end

function exact_phi(x::Float64, y::Float64, center=center_default())
    r2, theta = polar_features(x, y, center)
    return (r2^2) * cos(3 * theta)
end

function diffusivity_problem2(x::Float64, y::Float64, center=center_default())
    r2, _ = polar_features(x, y, center)
    return 1.0 - r2
end

function forcing_problem2(x::Float64, y::Float64, center=center_default())
    r2, theta = polar_features(x, y, center)
    return -(7.0 * r2 - 15.0 * r2^2) * cos(3 * theta)
end

function truncation_error_norms_l1_linf(solver, capacity, u_exact)
    centroids = capacity.C_ω
    u_exact_vals = [u_exact(c...) for c in centroids]
    residual = solver.A[1:end÷2, 1:end÷2] * u_exact_vals - solver.b[1:end÷2]
    cell_types = capacity.cell_types

    idx_all = findall((cell_types .== 1) .| (cell_types .== -1))
    idx_full = findall(cell_types .== 1)
    idx_cut = findall(cell_types .== -1)

    l1_all = lp_norm_subset(residual, idx_all, 1, capacity)
    l1_full = lp_norm_subset(residual, idx_full, 1, capacity)
    l1_cut = lp_norm_subset(residual, idx_cut, 1, capacity)

    linf_all = lp_norm_subset(residual, idx_all, Inf, capacity)
    linf_full = lp_norm_subset(residual, idx_full, Inf, capacity)
    linf_cut = lp_norm_subset(residual, idx_cut, Inf, capacity)

    return (
        l1_all = l1_all,
        l1_full = l1_full,
        l1_cut = l1_cut,
        linf_all = linf_all,
        linf_full = linf_full,
        linf_cut = linf_cut
    )
end

function run_truncation_convergence(
    nx_list::Vector{Int},
    ny_list::Vector{Int},
    center::Tuple{Float64,Float64},
    u_exact::Function;
    source::Function,
    diffusivity::Function,
    lx::Float64 = 1.0,
    ly::Float64 = 1.0,
    mean = nothing
)
    h_vals = Float64[]
    err_l1_all = Float64[]
    err_l1_full = Float64[]
    err_l1_cut = Float64[]
    err_l1_empty = Float64[]
    err_linf_all = Float64[]
    err_linf_full = Float64[]
    err_linf_cut = Float64[]
    err_linf_empty = Float64[]
    trunc_l1_all = Float64[]
    trunc_l1_full = Float64[]
    trunc_l1_cut = Float64[]
    trunc_linf_all = Float64[]
    trunc_linf_full = Float64[]
    trunc_linf_cut = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

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
        bc = Dirichlet((x,y,_=0)->u_exact(x,y))
        solver = if mean === nothing
            DiffusionSteadyMonoVariable(phase, bc_b, bc; mean=:harmonic)
        else
            DiffusionSteadyMonoVariable(phase, bc_b, bc; mean=:arithmetic)
        end
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err_l1, full_err_l1, cut_err_l1, empty_err_l1 =
            check_convergence(u_exact, solver, capacity, 1, false)
        _, _, global_err_linf, full_err_linf, cut_err_linf, empty_err_linf =
            check_convergence(u_exact, solver, capacity, Inf, false)
        trunc = truncation_error_norms_l1_linf(solver, capacity, u_exact)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_l1_all, global_err_l1)
        push!(err_l1_full, full_err_l1)
        push!(err_l1_cut, cut_err_l1)
        push!(err_l1_empty, empty_err_l1)
        push!(err_linf_all, global_err_linf)
        push!(err_linf_full, full_err_linf)
        push!(err_linf_cut, cut_err_linf)
        push!(err_linf_empty, empty_err_linf)
        push!(trunc_l1_all, trunc.l1_all)
        push!(trunc_l1_full, trunc.l1_full)
        push!(trunc_l1_cut, trunc.l1_cut)
        push!(trunc_linf_all, trunc.linf_all)
        push!(trunc_linf_full, trunc.linf_full)
        push!(trunc_linf_cut, trunc.linf_cut)
        push!(inside_cells, count_inside_cells(capacity))
        dx = lx / nx
        dy = ly / ny
        coverage_x = ceil(Int, 2 * star_max_radius() / dx)
        coverage_y = ceil(Int, 2 * star_max_radius() / dy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
    end

    return (
        h_vals = h_vals,
        err_l1_all = err_l1_all,
        err_l1_full = err_l1_full,
        err_l1_cut = err_l1_cut,
        err_l1_empty = err_l1_empty,
        err_linf_all = err_linf_all,
        err_linf_full = err_linf_full,
        err_linf_cut = err_linf_cut,
        err_linf_empty = err_linf_empty,
        trunc_l1_all = trunc_l1_all,
        trunc_l1_full = trunc_l1_full,
        trunc_l1_cut = trunc_l1_cut,
        trunc_linf_all = trunc_linf_all,
        trunc_linf_full = trunc_linf_full,
        trunc_linf_cut = trunc_linf_cut,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim
    )
end

function make_truncation_dataframe(method_name, data)
    pair_err_l1_all = compute_pairwise_orders(data.h_vals, data.err_l1_all)
    pair_err_l1_full = compute_pairwise_orders(data.h_vals, data.err_l1_full)
    pair_err_l1_cut = compute_pairwise_orders(data.h_vals, data.err_l1_cut)
    pair_err_linf_all = compute_pairwise_orders(data.h_vals, data.err_linf_all)
    pair_err_linf_full = compute_pairwise_orders(data.h_vals, data.err_linf_full)
    pair_err_linf_cut = compute_pairwise_orders(data.h_vals, data.err_linf_cut)

    pair_l1_all = compute_pairwise_orders(data.h_vals, data.trunc_l1_all)
    pair_l1_full = compute_pairwise_orders(data.h_vals, data.trunc_l1_full)
    pair_l1_cut = compute_pairwise_orders(data.h_vals, data.trunc_l1_cut)
    pair_linf_all = compute_pairwise_orders(data.h_vals, data.trunc_linf_all)
    pair_linf_full = compute_pairwise_orders(data.h_vals, data.trunc_linf_full)
    pair_linf_cut = compute_pairwise_orders(data.h_vals, data.trunc_linf_cut)

    return DataFrame(
        method = fill(method_name, length(data.h_vals)),
        h = data.h_vals,
        inside_cells = data.inside_cells,
        inside_cells_by_dim = data.inside_cells_by_dim,
        err_l1_all = data.err_l1_all,
        err_l1_full = data.err_l1_full,
        err_l1_cut = data.err_l1_cut,
        err_l1_empty = data.err_l1_empty,
        err_linf_all = data.err_linf_all,
        err_linf_full = data.err_linf_full,
        err_linf_cut = data.err_linf_cut,
        err_linf_empty = data.err_linf_empty,
        pair_order_err_l1_all = pair_err_l1_all,
        pair_order_err_l1_full = pair_err_l1_full,
        pair_order_err_l1_cut = pair_err_l1_cut,
        pair_order_err_linf_all = pair_err_linf_all,
        pair_order_err_linf_full = pair_err_linf_full,
        pair_order_err_linf_cut = pair_err_linf_cut,
        trunc_l1_all = data.trunc_l1_all,
        trunc_l1_full = data.trunc_l1_full,
        trunc_l1_cut = data.trunc_l1_cut,
        trunc_linf_all = data.trunc_linf_all,
        trunc_linf_full = data.trunc_linf_full,
        trunc_linf_cut = data.trunc_linf_cut,
        pair_order_l1_all = pair_l1_all,
        pair_order_l1_full = pair_l1_full,
        pair_order_l1_cut = pair_l1_cut,
        pair_order_linf_all = pair_linf_all,
        pair_order_linf_full = pair_linf_full,
        pair_order_linf_cut = pair_linf_cut
    )
end

function write_truncation_csv(method_name_default, data_default, method_name_arith, data_arith; csv_path=nothing)
    df_default = make_truncation_dataframe(method_name_default, data_default)
    df_arith = make_truncation_dataframe(method_name_arith, data_arith)
    df = vcat(df_default, df_arith)

    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "scalar") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "JohansenColella_P2_Truncation_Means.csv") : csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing, ny_list=nothing)
    center = center_default()
    nx_vals = isnothing(nx_list) ? [40, 80, 160, 320, 640] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list
    u_exact = (x,y) -> exact_phi(x,y,center)

    data_default = run_truncation_convergence(
        nx_vals, ny_vals, center, u_exact;
        source = (x,y,_)->forcing_problem2(x,y,center),
        diffusivity = (x,y,_)->diffusivity_problem2(x,y,center),
        lx = 1.0,
        ly = 1.0
    )

    data_arith = run_truncation_convergence(
        nx_vals, ny_vals, center, u_exact;
        source = (x,y,_)->forcing_problem2(x,y,center),
        diffusivity = (x,y,_)->diffusivity_problem2(x,y,center),
        lx = 1.0,
        ly = 1.0,
        mean = :arithmetic
    )

    csv_info = write_truncation_csv("default", data_default, "arithmetic", data_arith; csv_path=csv_path)
    return (
        data_default = data_default,
        data_arithmetic = data_arith,
        csv_path = csv_info.csv_path,
        table = csv_info.table
    )
end

results = main()

@testset "Johansen-Colella Problem 2 truncation mean comparison" begin
    @test length(results.data_default.h_vals) == length(results.data_arithmetic.h_vals)
    @test results.data_default.h_vals[1] > results.data_default.h_vals[end]
    @test isfile(results.csv_path)
end
