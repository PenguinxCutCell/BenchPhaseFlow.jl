using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
Liu–Fedkiw diphasic Poisson benchmark in 2D with a circular interface.

Domain: [0, 1] x [0, 1]. Interface: circle centered at (0.5, 0.5) with radius
0.25. Interior exact solution u1(x, y) = exp(-x^2 - y^2); exterior u2(x, y) = 0.
Diffusivity: D1 = 2 (interior), D2 = 1 (exterior).

Source term:
    f1(x, y) = 8(x^2 + y^2 - 1)exp(-x^2 - y^2)  (interior)
    f2(x, y) = 0                                (exterior)

Interface jumps (phase 2 minus phase 1):
    [u] = -exp(-x^2 - y^2)
    [D∇u·n] = 8(2x^2 + 2y^2 - x - y)exp(-x^2 - y^2)
Outer boundaries: Dirichlet 0.
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

center_default() = (0.5, 0.5)
radius_default() = 0.25
domain_lengths() = (1.0, 1.0)
domain_origin() = (0.0, 0.0)

u_interior(x, y) = exp(-(x^2 + y^2))
u_exterior(x, y) = 0.0

f_interior(x, y, _=0.0) = 8.0 * ((x^2 + y^2) - 1.0) * exp(-(x^2 + y^2))
f_exterior(x, y, _=0.0) = 0.0

D_interior(x, y, _=0.0) = 2.0
D_exterior(x, y, _=0.0) = 1.0

jump_scalar(x, y, _=0.0) = -exp(-(x^2 + y^2))
jump_flux(x, y, _=0.0) = 8.0 * (2.0 * x^2 + 2.0 * y^2 - x - y) * exp(-(x^2 + y^2))

function build_poisson_circle_solver(mesh; center=center_default(), radius=radius_default())
    body_out = (x, y, _=0.0) -> hypot(x - center[1], y - center[2]) - radius
    body_in  = (x, y, _=0.0) -> radius - hypot(x - center[1], y - center[2])

    capacity1 = Capacity(body_in, mesh)   # interior
    capacity2 = Capacity(body_out, mesh)  # exterior
    operator1 = DiffusionOps(capacity1)
    operator2 = DiffusionOps(capacity2)

    bc_b = BorderConditions(Dict(
        :left => Dirichlet(0.0),
        :right => Dirichlet(0.0),
        :top => Dirichlet(0.0),
        :bottom => Dirichlet(0.0)
    ))

    ic = InterfaceConditions(
        ScalarJump(1.0, 1.0, jump_scalar),
        FluxJump(2.0, 1.0, jump_flux)
    )

    phase1 = Phase(capacity1, operator1, f_interior, D_interior)
    phase2 = Phase(capacity2, operator2, f_exterior, D_exterior)

    solver = DiffusionSteadyDiph(phase1, phase2, bc_b, ic)
    return solver, capacity1, capacity2
end

function run_liufedkiw_poisson2d(
    nx_list::Vector{Int};
    lx::Float64 = domain_lengths()[1],
    ly::Float64 = domain_lengths()[2],
    center::Tuple{Float64,Float64} = center_default(),
    radius::Float64 = radius_default(),
    norm::Real = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    phase1_all_errs = Float64[]
    phase1_full_errs = Float64[]
    phase1_cut_errs = Float64[]
    phase1_empty_errs = Float64[]
    phase2_all_errs = Float64[]
    phase2_full_errs = Float64[]
    phase2_cut_errs = Float64[]
    phase2_empty_errs = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()
    inside_cells_phase1 = Int[]
    inside_cells_phase2 = Int[]

    for nx in nx_list
        ny = nx
        mesh = Penguin.Mesh((nx, ny), (lx, ly), domain_origin())
        solver, capacity1, capacity2 = build_poisson_circle_solver(mesh; center=center, radius=radius)

        solve_DiffusionSteadyDiph!(solver; method=Base.:\)
        push!(solver.states, solver.x)

        _, _, global_errs, full_errs, cut_errs, empty_errs =
            check_convergence_diph(u_interior, u_exterior, solver, capacity1, capacity2, norm, relative)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, global_errs[3])
        push!(err_full_vals, full_errs[3])
        push!(err_cut_vals, cut_errs[3])
        push!(err_empty_vals, empty_errs[3])
        push!(phase1_all_errs, global_errs[1])
        push!(phase2_all_errs, global_errs[2])
        push!(phase1_full_errs, full_errs[1])
        push!(phase2_full_errs, full_errs[2])
        push!(phase1_cut_errs, cut_errs[1])
        push!(phase2_cut_errs, cut_errs[2])
        push!(phase1_empty_errs, empty_errs[1])
        push!(phase2_empty_errs, empty_errs[2])

        inside1 = count_inside_cells(capacity1)
        inside2 = count_inside_cells(capacity2)
        push!(inside_cells, inside1 + inside2)
        push!(inside_cells_by_dim, [inside1, inside2])
        push!(inside_cells_phase1, inside1)
        push!(inside_cells_phase2, inside2)
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full_vals,
        err_cut_vals = err_cut_vals,
        err_empty_vals = err_empty_vals,
        phase1_all_errs = phase1_all_errs,
        phase1_full_errs = phase1_full_errs,
        phase1_cut_errs = phase1_cut_errs,
        phase1_empty_errs = phase1_empty_errs,
        phase2_all_errs = phase2_all_errs,
        phase2_full_errs = phase2_full_errs,
        phase2_cut_errs = phase2_cut_errs,
        phase2_empty_errs = phase2_empty_errs,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        inside_cells_phase1 = inside_cells_phase1,
        inside_cells_phase2 = inside_cells_phase2,
        orders = compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals),
        norm = norm
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_diphasic_convergence_dataframe(method_name, data)
    results_dir = isnothing(csv_path) ?
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic", "LiuFedkiw") :
        dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ?
        joinpath(results_dir, "$(method_name)_Convergence.csv") :
        csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing)
    nx_vals = isnothing(nx_list) ? [16, 32, 64, 128] : nx_list
    data = run_liufedkiw_poisson2d(nx_vals)
    csv_info = write_convergence_csv("LiuFedkiw_Poisson2D_Circle", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Liu-Fedkiw Diphasic Poisson 2D Circle" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end
