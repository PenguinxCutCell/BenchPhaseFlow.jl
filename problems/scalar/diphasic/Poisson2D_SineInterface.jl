using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test
using CairoMakie

"""
Diphasic 2D Poisson convergence test with a sine-wave interface.

Domain: [0, 1] x [0, 1]. Interface y = 0.5 + 0.1 * sin(2πx).
Phase 1: y <= interface; Phase 2: y >= interface.
Diffusivities: D1 = 1, D2 = 2.

Manufactured solution (zero on the box boundary):
    φ(x, y) = sin(πx) * sin(πy)
    u1 = a * φ, u2 = b * φ with b = a * D1 / D2
This enforces flux continuity (β2 ∂n u2 = β1 ∂n u1) while allowing a scalar jump
[u] = (b - a) φ. Set a = 1 here.

Laplacian: Δφ = -2π^2 φ, so forcing per phase:
    f1 = 2π^2 * D1 * a * φ
    f2 = 2π^2 * D2 * b * φ = 2π^2 * D1 * a * φ (same value).

Interface conditions:
    ScalarJump = (b - a) * φ
    FluxJump   = 0
Boundary: Dirichlet exact solution (φ = 0 on all borders).
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

interface_y(x) = 0.51 + 0.1 * sin(2π * x)
domain_origin() = (0.0, 0.0)
domain_lengths() = (1.0, 1.0)

φ(x, y) = sin(π * x) * sin(π * y)
const A_SOL = 1.0
const D1_CONST = 1.0
const D2_CONST = 1.0
const B_SOL = A_SOL * D1_CONST / D2_CONST

u1_exact(x, y) = A_SOL * φ(x, y)
u2_exact(x, y) = B_SOL * φ(x, y)
f_rhs_common(x, y, _=0.0) = 2π^2 * D1_CONST * A_SOL * φ(x, y) # same for both
D1(x, y, _=0.0) = D1_CONST
D2(x, y, _=0.0) = D2_CONST

function build_sine_interface_capacities(mesh)
    body1 = (x, y, _=0.0) -> interface_y(x) - y   # inside (phase 1) below curve
    body2 = (x, y, _=0.0) -> y - interface_y(x)   # outside (phase 2) above curve
    capacity1 = Capacity(body1, mesh; method="VOFI", integration_method=:vofijul)
    capacity2 = Capacity(body2, mesh; method="VOFI", integration_method=:vofijul)
    operator1 = DiffusionOps(capacity1)
    operator2 = DiffusionOps(capacity2)
    return capacity1, capacity2, operator1, operator2
end

function run_poisson_sine_interface(
    nx_list::Vector{Int};
    lx::Float64 = domain_lengths()[1],
    ly::Float64 = domain_lengths()[2],
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
        capacity1, capacity2, operator1, operator2 = build_sine_interface_capacities(mesh)

        bc_func = (x, y, _=0.0) -> u1_exact(x, y)
        bc_b = BorderConditions(Dict(
            :left => Dirichlet(0.0),
            :right => Dirichlet(0.0),
            :top => Dirichlet(0.0),
            :bottom => Dirichlet(0.0)
        ))

        jump_val = (x, y, _=0.0) -> (B_SOL - A_SOL) * φ(x, y)
        ic = InterfaceConditions(
            ScalarJump(1.0, 1.0, 0.0),
            FluxJump(D1_CONST, D2_CONST, 0.0)
        )

        phase1 = Phase(capacity1, operator1, f_rhs_common, D1)
        phase2 = Phase(capacity2, operator2, f_rhs_common, D2)

        solver = DiffusionSteadyDiph(phase1, phase2, bc_b, ic)
        solve_DiffusionSteadyDiph!(solver; method=Base.:\)
        push!(solver.states, solver.x)

        
        (u1_ana, u2_ana),(u1_num, u2_num), global_errs, full_errs, cut_errs, empty_errs =
            check_convergence_diphh(u1_exact, u2_exact, solver, capacity1, capacity2, norm, relative)

        # Compute errors
        err1 = u1_ana .- u1_num
        err2 = u2_ana .- u2_num
        
        # Get cell types for each phase
        cell_types1 = capacity1.cell_types
        cell_types2 = capacity2.cell_types
        
        # Filter indices by cell type for phase 1
        idx_all1 = findall((cell_types1 .== 1) .| (cell_types1 .== -1))
        idx_full1 = findall(cell_types1 .== 1)
        idx_cut1 = findall(cell_types1 .== -1)
        idx_empty1 = findall(cell_types1 .== 0)

        # Filter indices by cell type for phase 2
        idx_all2 = findall((cell_types2 .== 1) .| (cell_types2 .== -1))
        idx_full2 = findall(cell_types2 .== 1)
        idx_cut2 = findall(cell_types2 .== -1)
        idx_empty2 = findall(cell_types2 .== 0) 

        #Extract only full cells for debugging plots
        err1_full = zeros(size(err1))
        err2_full = zeros(size(err2))
        for i in idx_all1
            err1_full[i] = err1[i]
        end
        for i in idx_full2
            err2_full[i] = err2[i]
        end
        

        #Debugging plots
        fig = Figure(resolution = (800, 400))
        ax1 = Axis(fig[1, 1], title = "Phase 1 Numerical", xlabel = "x", ylabel = "y")
        hm = heatmap!(ax1, mesh.nodes[1], mesh.nodes[2], reshape(u1_num, nx+1, ny+1)')
        Colorbar(fig[1, 2], hm; label = "u1_num")
        ax2 = Axis(fig[2, 1], title = "Phase 1 Analytical", xlabel = "x", ylabel = "y")
        hm2 = heatmap!(ax2, mesh.nodes[1], mesh.nodes[2], reshape(u1_ana, nx+1, ny+1)')
        Colorbar(fig[2, 2], hm2; label = "u1_ana")
        ax3 = Axis(fig[3, 1], title = "Phase 1 Error", xlabel = "x", ylabel = "y")
        hm3 = heatmap!(ax3, mesh.nodes[1], mesh.nodes[2], reshape(err1_full, nx+1, ny+1)')
        Colorbar(fig[3, 2], hm3; label = "Error")
        display(fig)

        readline()

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
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic") :
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
    data = run_poisson_sine_interface(nx_vals)
    csv_info = write_convergence_csv("Diph_Poisson2D_SineInterface", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Diphasic Poisson 2D sine interface" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end
