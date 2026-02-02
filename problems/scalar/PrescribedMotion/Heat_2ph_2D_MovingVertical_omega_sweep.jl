using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using Test

"""
Diphasic 2D heat equation with a vertically oscillating interface.
Sweeps ω from 2π to 32π and writes convergence data to CSV for each case.
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

const DEFAULT_OMEGA_SWEEP = [4π, 8π, 16π, 32π]

struct MovingVerticalParams
    lx::Float64
    ly::Float64
    x0::Float64
    y0::Float64
    Tend::Float64
    D_plus::Float64
    D_minus::Float64
    cp_plus::Float64
    cp_minus::Float64
    s0::Float64
    A::Float64
    ω::Float64
end

MovingVerticalParams(; lx=4.0, ly=4.0, x0=0.0, y0=0.0, Tend=0.1,
                    D_plus=1.0, D_minus=1.0, cp_plus=1.0, cp_minus=1.0,
                    s0=2.0, A=0.1, ω=32π) = MovingVerticalParams(
    lx, ly, x0, y0, Tend, D_plus, D_minus, cp_plus, cp_minus, s0, A, ω
)

s_val(t, p::MovingVerticalParams) = p.s0 + p.A * sin(p.ω * t)
sdot_val(t, p::MovingVerticalParams) = p.A * p.ω * cos(p.ω * t)

g1(x, p::MovingVerticalParams) = (x - p.x0) * (p.lx - (x - p.x0))
g2(y, p::MovingVerticalParams) = (y - p.y0) * (p.ly - (y - p.y0))
g(x, y, p::MovingVerticalParams) = g1(x, p) * g2(y, p)
θ(t) = exp(-t)

function u1_exact(p::MovingVerticalParams)
    return (x, y, t) -> x >= s_val(t, p) ? p.D_minus * (x - s_val(t, p)) * g(x, y, p) * θ(t) : 0.0
end

function u2_exact(p::MovingVerticalParams)
    return (x, y, t) -> x < s_val(t, p) ? p.D_plus * (x - s_val(t, p)) * g(x, y, p) * θ(t) : 0.0
end

function source_terms(p::MovingVerticalParams)
    f1 = (x, y, z, t) -> begin
        if x < s_val(t, p)
            return 0.0
        end
        time_term = -p.cp_plus * p.D_minus * θ(t) * g(x, y, p) * (cos(p.ω * t) + (x - s_val(t, p)))
        lap_g = (8.0 - 6.0 * x + 2.0 * s_val(t, p)) * g2(y, p) - 2.0 * (x - s_val(t, p)) * g1(x, p)
        space_term = -p.D_plus * θ(t) * lap_g
        return time_term + space_term
    end

    f2 = (x, y, z, t) -> begin
        if x >= s_val(t, p)
            return 0.0
        end
        time_term = -p.cp_minus * p.D_plus * θ(t) * g(x, y, p) * (cos(p.ω * t) + (x - s_val(t, p)))
        lap_g = (8.0 - 6.0 * x + 2.0 * s_val(t, p)) * g2(y, p) - 2.0 * (x - s_val(t, p)) * g1(x, p)
        space_term = -p.D_minus * θ(t) * lap_g
        return time_term + space_term
    end

    return f1, f2
end

function run_moving_vertical_convergence(
    nx_list::Vector{Int};
    params::MovingVerticalParams=MovingVerticalParams(),
    norm::Real=2,
    relative::Bool=false
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

    u1 = u1_exact(params)
    u2 = u2_exact(params)
    f1, f2 = source_terms(params)

    for nx in nx_list
        ny = nx
        mesh = Penguin.Mesh((nx, ny), (params.lx, params.ly), (params.x0, params.y0))
        Δt = 0.5 * (params.lx / nx)^2
        STmesh = Penguin.SpaceTimeMesh(mesh, [0.0, Δt], tag=mesh.tag)

        body = (x, y, t) -> x - s_val(t, params)
        body_c = (x, y, t) -> s_val(t, params) - x

        capacity = Capacity(body, STmesh)
        capacity_c = Capacity(body_c, STmesh)
        operator = DiffusionOps(capacity)
        operator_c = DiffusionOps(capacity_c)

        bc = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left => bc,
            :right => bc,
            :top => bc,
            :bottom => bc
        ))

        ic = InterfaceConditions(
            ScalarJump(1.0, 1.0, 0.0),
            FluxJump(params.D_plus, params.D_minus, 0.0)
        )

        K1_func = (x, y, z) -> params.D_plus
        K2_func = (x, y, z) -> params.D_minus

        phase1 = Phase(capacity, operator, f1, K1_func)
        phase2 = Phase(capacity_c, operator_c, f2, K2_func)

        ndofs = (nx + 1) * (ny + 1)
        u0ₒ1 = zeros(ndofs)
        u0ᵧ1 = zeros(ndofs)
        u0ₒ2 = zeros(ndofs)
        u0ᵧ2 = zeros(ndofs)

        for j in 1:ny+1, i in 1:nx+1
            idx = (j - 1) * (nx + 1) + i
            x = params.x0 + (i - 1) * params.lx / nx
            y = params.y0 + (j - 1) * params.ly / ny
            if x >= s_val(0.0, params)
                u0ₒ1[idx] = u1(x, y, 0.0)
                u0ᵧ1[idx] = u1(x, y, 0.0)
            else
                u0ₒ2[idx] = u2(x, y, 0.0)
                u0ᵧ2[idx] = u2(x, y, 0.0)
            end
        end

        u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

        solver = MovingDiffusionUnsteadyDiph(phase1, phase2, bc_b, ic, Δt, u0, mesh, "BE")
        solve_MovingDiffusionUnsteadyDiph!(solver, phase1, phase2, body, body_c, Δt, params.Tend, bc_b, ic, mesh, "BE"; method=Base.:\)

        body_tend = (x, y, _=0) -> x - s_val(params.Tend, params)
        body_tend_c = (x, y, _=0) -> s_val(params.Tend, params) - x
        capacity_tend = Capacity(body_tend, mesh; compute_centroids=false)
        capacity_tend_c = Capacity(body_tend_c, mesh; compute_centroids=false)

        _, _, global_errs, full_errs, cut_errs, empty_errs =
            check_convergence_diph((x, y)->u1(x, y, params.Tend), (x, y)->u2(x, y, params.Tend),
                                   solver, capacity_tend, capacity_tend_c, norm, relative)

        push!(h_vals, params.lx / nx)
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

        inside1 = count_inside_cells(capacity_tend)
        inside2 = count_inside_cells(capacity_tend_c)
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
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic", "prescribed_motion") :
        dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ?
        joinpath(results_dir, "$(method_name)_Convergence.csv") :
        csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function _rebuild_params(params::MovingVerticalParams; ω=params.ω)
    return MovingVerticalParams(
        params.lx, params.ly, params.x0, params.y0, params.Tend,
        params.D_plus, params.D_minus, params.cp_plus, params.cp_minus,
        params.s0, params.A, ω
    )
end

function main(;
    csv_path=nothing,
    csv_dir=csv_path,
    nx_list=nothing,
    omega_list=nothing,
    params::MovingVerticalParams=MovingVerticalParams()
)
    nx_vals = isnothing(nx_list) ? [4, 8, 16, 32, 64, 128] : nx_list
    omegas = isnothing(omega_list) ? DEFAULT_OMEGA_SWEEP : omega_list

    if length(omegas) > 1 && !isnothing(csv_dir) && endswith(lowercase(csv_dir), ".csv")
        error("Provide a directory path (or nothing) for csv_dir when running multiple cases.")
    end

    results = NamedTuple[]

    for ω in omegas
        params_case = _rebuild_params(params; ω=ω)
        data = run_moving_vertical_convergence(nx_vals; params=params_case)
        method_name = "Heat_2ph_2D_MovingVertical_omega$(params_case.ω)"
        csv_info = write_convergence_csv(method_name, data; csv_path=csv_dir)
        push!(results, (
            ω = params_case.ω,
            params = params_case,
            data = data,
            csv_path = csv_info.csv_path,
            table = csv_info.table
        ))
    end

    return results
end

results = main()

@testset "Diphasic moving vertical interface convergence (omega sweep)" begin
    @test length(results) == length(DEFAULT_OMEGA_SWEEP)
    for res in results
        orders = res.data.orders
        @test !isnan(orders.all)
        @test length(res.data.h_vals) == length(res.data.err_vals)
        @test res.data.h_vals[1] > res.data.h_vals[end]
        @test minimum(res.data.err_vals) < maximum(res.data.err_vals)
        @test isfile(res.csv_path)
    end
end
