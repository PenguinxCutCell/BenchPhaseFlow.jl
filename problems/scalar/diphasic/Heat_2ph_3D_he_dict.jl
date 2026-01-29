using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using QuadGK
using CSV
using Test

"""
Diphasic 3D diffusion benchmark (spherical interface) with analytical radial solution.

This variant sweeps Henry coefficients and runs the mesh convergence study for
each value. The analytical reference includes the Henry jump and a simple
interface consistency check is reported per case.
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

const DEFAULT_HENRY_SWEEP = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]

struct Heat2Ph3DParams
    lx::Float64
    ly::Float64
    lz::Float64
    center::Tuple{Float64,Float64,Float64}
    radius::Float64
    Tend::Float64
    Dg::Float64
    Dl::Float64
    He::Float64
    c0::Float64
    X_MAX::Float64
end

Heat2Ph3DParams(; lx=8.0, ly=8.0, lz=8.0, center=(2.0, 2.0, 2.0),
                radius=1.0, Tend=0.1, Dg=1.0, Dl=1.0, He=1.0,
                c0=1.0, X_MAX=1000.0) =
    Heat2Ph3DParams(lx, ly, lz, center, radius, Tend, Dg, Dl, He, c0, X_MAX)

function make_exact_solution(params::Heat2Ph3DParams)
    σ = sqrt(params.Dl / params.Dg)
    Q = (params.Dl / params.Dg) * σ
    L = (params.Dl - params.Dg) / params.Dg

    t1_cache = Dict{Tuple{Float64,Float64},Float64}()
    t2_cache = Dict{Tuple{Float64,Float64},Float64}()

    denom(u) = (u * cos(u) + L * sin(u))^2 + (Q * u * sin(u))^2

    function integrand_core(u, r, t)
        num = (sin(u) - u * cos(u)) * sin(u * r / params.radius)
        expo = exp(-params.Dg * u^2 * t / params.radius^2)
        return num * expo / denom(u)
    end

    function F_shell(u, r)
        arg = u * (r - params.radius) / (σ * params.radius)
        return (u * cos(u) + L * sin(u)) * sin(arg) + Q * u * sin(u) * cos(arg)
    end

    function integrand_shell(u, r, t)
        num = (sin(u) - u * cos(u)) * F_shell(u, r)
        expo = exp(-params.Dg * u^2 * t / params.radius^2)
        u_safe = max(u, 1e-12)
        return num * expo * σ / (u_safe * denom(u))
    end

    function T_core_exact(t, r)
        r_eff = max(r, 1e-8)
        key = (round(t, sigdigits=6), round(r_eff, sigdigits=6))
        return get!(t1_cache, key) do
            val, _ = quadgk(u -> integrand_core(u, r_eff, t), 0.0, params.X_MAX;
                            atol=1e-6, rtol=1e-6, maxevals=10_000)
            (2 * Q * params.He / (pi * r_eff)) * val * params.c0
        end
    end

    function T_shell_exact(t, r)
        r_eff = max(r, 1e-8)
        key = (round(t, sigdigits=6), round(r_eff, sigdigits=6))
        return get!(t2_cache, key) do
            val, _ = quadgk(u -> integrand_shell(u, r_eff, t), 0.0, params.X_MAX;
                            atol=1e-6, rtol=1e-6, maxevals=10_000)
            (2 * params.He / (pi * r_eff)) * val * params.c0
        end
    end

    exact_field = function (t, x, y, z)
        r = sqrt((x - params.center[1])^2 + (y - params.center[2])^2 + (z - params.center[3])^2)
        r < params.radius ? T_core_exact(t, r) : T_shell_exact(t, r)
    end

    return (T_core_exact=T_core_exact, T_shell_exact=T_shell_exact, exact_field=exact_field)
end

f_zero(x, y, z, t) = 0.0

function build_bubble_components(mesh, params::Heat2Ph3DParams)
    body_out = (x, y, z, _=0.0) -> sqrt((x - params.center[1])^2 +
                                        (y - params.center[2])^2 +
                                        (z - params.center[3])^2) - params.radius
    body_in  = (x, y, z, _=0.0) -> params.radius - sqrt((x - params.center[1])^2 +
                                                        (y - params.center[2])^2 +
                                                        (z - params.center[3])^2)

    capacity1 = Capacity(body_in, mesh; method="VOFI")
    capacity2 = Capacity(body_out, mesh; method="VOFI")
    operator1 = DiffusionOps(capacity1)
    operator2 = DiffusionOps(capacity2)

    bc_b = BorderConditions(Dict{Symbol,AbstractBoundary}())

    ic = InterfaceConditions(
        ScalarJump(1.0, 1.0 / params.He, 0.0),
        FluxJump(params.Dg, params.Dl, 0.0)
    )

    phase1 = Phase(capacity1, operator1, f_zero, (x, y, z) -> params.Dg)
    phase2 = Phase(capacity2, operator2, f_zero, (x, y, z) -> params.Dl)

    return capacity1, capacity2, operator1, operator2, phase1, phase2, bc_b, ic
end

function interface_consistency(exact, params::Heat2Ph3DParams; t=params.Tend)
    core_val = exact.T_core_exact(t, params.radius)
    shell_val = exact.T_shell_exact(t, params.radius)
    mismatch = shell_val - params.He * core_val
    return (core=core_val, shell=shell_val, mismatch=mismatch)
end

function run_bubble_3d(
    nx_list::Vector{Int};
    params::Heat2Ph3DParams=Heat2Ph3DParams(),
    norm::Real = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    dt_vals = Float64[]
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

    exact = make_exact_solution(params)
    iface_check = interface_consistency(exact, params)

    for nx in nx_list
        ny = nx; nz = nx
        mesh = Penguin.Mesh((nx, ny, nz), (params.lx, params.ly, params.lz), (0.0, 0.0, 0.0))

        capacity1, capacity2, operator1, operator2, phase1, phase2, bc_b, ic =
            build_bubble_components(mesh, params)
        ndofs1 = length(capacity1.C_ω)
        ndofs2 = length(capacity2.C_ω)
        @assert ndofs1 == ndofs2 "Expected matching DOF counts for both phases"

        u0_phase1 = [exact.exact_field(0.0, c[1], c[2], c[3]) for c in capacity1.C_ω]
        u0_phase2 = [exact.exact_field(0.0, c[1], c[2], c[3]) for c in capacity2.C_ω]
        u0 = vcat(u0_phase1, u0_phase1, u0_phase2, u0_phase2)

        Δt = 0.5 * (minimum((params.lx / nx, params.ly / ny, params.lz / nz))^2) / max(params.Dg, params.Dl)
        push!(dt_vals, Δt)

        solver = DiffusionUnsteadyDiph(phase1, phase2, bc_b, ic, Δt, u0, "BE")
        solve_DiffusionUnsteadyDiph!(solver, phase1, phase2, Δt, params.Tend, bc_b, ic, "BE"; method=Base.:\)
        push!(solver.states, solver.x)

        u_exact = (x, y, z) -> exact.exact_field(params.Tend, x, y, z)
        _, _, global_errs, full_errs, cut_errs, empty_errs =
            check_convergence_diphh(u_exact, u_exact, solver, capacity1, capacity2, norm, relative)

        push!(h_vals, minimum((params.lx / nx, params.ly / ny, params.lz / nz)))
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
        dt_vals = dt_vals,
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
        interface_core = iface_check.core,
        interface_shell = iface_check.shell,
        interface_mismatch = iface_check.mismatch,
        norm = norm
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_diphasic_convergence_dataframe(method_name, data)
    df.interface_core = fill(data.interface_core, nrow(df))
    df.interface_shell = fill(data.interface_shell, nrow(df))
    df.interface_mismatch = fill(data.interface_mismatch, nrow(df))
    results_dir = isnothing(csv_path) ?
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic", "he_sweep_3d") :
        dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ?
        joinpath(results_dir, "$(method_name)_Convergence.csv") :
        csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function _rebuild_params(params::Heat2Ph3DParams; He=params.He)
    return Heat2Ph3DParams(
        params.lx, params.ly, params.lz, params.center, params.radius,
        params.Tend, params.Dg, params.Dl, He, params.c0, params.X_MAX
    )
end

function main(;
    csv_path=nothing,
    csv_dir=csv_path,
    nx_list=nothing,
    He_list=nothing,
    params::Heat2Ph3DParams=Heat2Ph3DParams()
)
    nx_vals = isnothing(nx_list) ? [4, 8, 16, 31, 42, 53] : nx_list
    He_vals = isnothing(He_list) ? DEFAULT_HENRY_SWEEP : He_list

    if length(He_vals) > 1 && !isnothing(csv_dir) && endswith(lowercase(csv_dir), ".csv")
        error("Provide a directory path (or nothing) for csv_dir when running multiple cases.")
    end

    results = NamedTuple[]

    for He in He_vals
        params_case = _rebuild_params(params; He=He)
        data = run_bubble_3d(nx_vals; params=params_case)
        method_name = "Diph_3D_Bubble_He$(params_case.He)"
        csv_info = write_convergence_csv(method_name, data; csv_path=csv_dir)
        push!(results, (
            He = params_case.He,
            params = params_case,
            data = data,
            csv_path = csv_info.csv_path,
            table = csv_info.table
        ))
    end

    return results
end

results = main()

@testset "Diphasic 3D bubble convergence (He sweep)" begin
    @test length(results) == length(DEFAULT_HENRY_SWEEP)
    for res in results
        orders = res.data.orders
        @test !isnan(orders.all)
        @test length(res.data.h_vals) == length(res.data.err_vals)
        @test res.data.h_vals[1] > res.data.h_vals[end]
        @test minimum(res.data.err_vals) < maximum(res.data.err_vals)
        @test isfinite(res.data.interface_mismatch)
        @test isfile(res.csv_path)
    end
end
