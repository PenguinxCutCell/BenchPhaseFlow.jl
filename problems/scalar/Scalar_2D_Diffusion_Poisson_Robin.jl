using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
Steady 2D Poisson problem with an embedded circular Robin boundary.

The domain is [0, 4]^2 with a circular obstacle of radius `ly / 4` centered at
`(lx/2, ly/2) + (0.01, 0.01)`. Dirichlet zero data is imposed on the outer box
and Robin(1, 1, 1) data on the embedded boundary. The manufactured solution
`u(x, y) = 7/4 - ((x - cx)^2 + (y - cy)^2)/4` is used to measure both L2 and H1
(gradient) convergence.
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

domain_origin() = (0.0, 0.0)
domain_lengths() = (4.0, 4.0)
default_radius(lx, ly) = ly / 4
default_center(lx, ly) = (lx / 2 + 0.01, ly / 2 + 0.01)

circle_level_set(x, y, radius, center) = sqrt((x - center[1])^2 + (y - center[2])^2) - radius
u_exact_disk(x, y, center) = 7.0 / 4.0 - ((x - center[1])^2 + (y - center[2])^2) / 4.0
grad_exact_disk(center) = (
    (x, y) -> -(x - center[1]) / 2,
    (x, y) -> -(y - center[2]) / 2
)

source_constant(x, y, _=0.0) = 1.0
diffusivity_constant(x, y, _=0.0) = 1.0

function run_poisson_robin_convergence(
    nx_list::Vector{Int},
    ny_list::Vector{Int};
    lx::Float64 = domain_lengths()[1],
    ly::Float64 = domain_lengths()[2],
    center::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    radius::Union{Nothing,Float64} = nothing,
    norm::Int = 2,
    relative::Bool = false
)
    c = isnothing(center) ? default_center(lx, ly) : center
    r = isnothing(radius) ? default_radius(lx, ly) : radius
    u_exact = (x, y) -> u_exact_disk(x, y, c)
    grad_exact = grad_exact_disk(c)

    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    h1_err_vals = Float64[]
    h1_full_vals = Float64[]
    h1_cut_vals = Float64[]
    h1_empty_vals = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    for (nx, ny) in zip(nx_list, ny_list)
        mesh = Penguin.Mesh((nx, ny), (lx, ly), domain_origin())
        capacity = Capacity((x, y, _=0.0) -> circle_level_set(x, y, r, c), mesh)
        operator = DiffusionOps(capacity)

        bc_interface = Robin(1.0, 1.0, 1.0)
        bc_outer = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left => bc_outer,
            :right => bc_outer,
            :top => bc_outer,
            :bottom => bc_outer
        ))
        phase = Phase(capacity, operator, source_constant, diffusivity_constant)
        solver = DiffusionSteadyMono(phase, bc_b, bc_interface)

        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, all_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)

        _, _, h1_all_err, h1_full_err, h1_cut_err, h1_empty_err =
            check_h1_convergence(grad_exact, solver, capacity, operator; p=norm, relative=relative)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, all_err)
        push!(err_full_vals, full_err)
        push!(err_cut_vals, cut_err)
        push!(err_empty_vals, empty_err)
        push!(h1_err_vals, h1_all_err)
        push!(h1_full_vals, h1_full_err)
        push!(h1_cut_vals, h1_cut_err)
        push!(h1_empty_vals, h1_empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        dx = lx / nx
        dy = ly / ny
        coverage_x = ceil(Int, 2 * r / dx)
        coverage_y = ceil(Int, 2 * r / dy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full_vals,
        err_cut_vals = err_cut_vals,
        err_empty_vals = err_empty_vals,
        h1_err_vals = h1_err_vals,
        h1_full_vals = h1_full_vals,
        h1_cut_vals = h1_cut_vals,
        h1_empty_vals = h1_empty_vals,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        orders = compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals),
        h1_orders = compute_orders(h_vals, h1_err_vals, h1_full_vals, h1_cut_vals),
        norm = norm
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_convergence_dataframe(method_name, data)

    h1_pair_all = compute_pairwise_orders(data.h_vals, data.h1_err_vals)
    h1_pair_full = compute_pairwise_orders(data.h_vals, data.h1_full_vals)
    h1_pair_cut = compute_pairwise_orders(data.h_vals, data.h1_cut_vals)

    df.h1_all_err = data.h1_err_vals
    df.h1_full_err = data.h1_full_vals
    df.h1_cut_err = data.h1_cut_vals
    df.h1_empty_err = data.h1_empty_vals
    df.pair_order_h1_all = h1_pair_all
    df.pair_order_h1_full = h1_pair_full
    df.pair_order_h1_cut = h1_pair_cut

    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "scalar") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "$(method_name)_Convergence.csv") : csv_path

    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing, ny_list=nothing)
    lx, ly = domain_lengths()
    nx_vals = isnothing(nx_list) ? [16, 32, 64, 128, 256, 512] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list

    data = run_poisson_robin_convergence(
        nx_vals, ny_vals;
        lx=lx, ly=ly,
        norm=2,
        relative=false
    )

    csv_info = write_convergence_csv("Scalar_2D_Diffusion_Poisson_Robin", data; csv_path=csv_path)
    return (data=data, csv_path=csv_info.csv_path, table=csv_info.table)
end

results = main()

@testset "Poisson Robin 2D convergence" begin
    orders = results.data.orders
    h1_orders = results.data.h1_orders
    @test orders.all > 1.0
    @test h1_orders.all > 0.5
    @test length(results.data.h_vals) == length(results.data.err_vals) == length(results.data.h1_err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test minimum(results.data.h1_err_vals) < maximum(results.data.h1_err_vals)
    @test isfile(results.csv_path)
end
