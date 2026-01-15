using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using DataFrames
using Test

"""
Parametric study for the steep k-profile disk Poisson problem.
Sweeps eps (steepness), k_max/k_min (jump), and mean type (arithmetic/harmonic).
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

function k_profile(r; k_min=1.0, k_max=10.0, r0=0.25, eps=0.025)
    delta_k = k_max - k_min
    return k_min + 0.5 * delta_k * (1.0 + tanh((r - r0) / eps))
end

function k_profile_prime(r; k_min=1.0, k_max=10.0, r0=0.25, eps=0.025)
    delta_k = k_max - k_min
    z = (r - r0) / eps
    sech2 = 1.0 / cosh(z)^2
    return 0.5 * delta_k / eps * sech2
end

function poly_real_part(x, y, m)
    z = complex(x, y)
    return real(z^m)
end

function exact_solution(x, y, r, R, m, center)
    xr = x - center[1]
    yr = y - center[2]
    return (1.0 - (r / R)^2) * poly_real_part(xr, yr, m)
end

function forcing_term(x, y, r, R, m, center; k_min=1.0, k_max=10.0, r0=0.25, eps=0.025)
    k_val = k_profile(r; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
    k_prime = k_profile_prime(r; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
    if r == 0.0
        return 0.0
    end

    xr = x - center[1]
    yr = y - center[2]
    p_m = poly_real_part(xr, yr, m)
    a_prime = m * r^(m - 1) - (m + 2) * r^(m + 1) / (R^2)
    return (4.0 * (m + 1) / (R^2)) * k_val * p_m - (k_prime * a_prime) * (p_m / r^m)
end

function run_disk_convergence(
    nx_list::Vector{Int};
    lx::Float64 = 1.0,
    ly::Float64 = 1.0,
    center::Tuple{Float64,Float64} = (0.5, 0.5),
    R::Float64 = 0.5,
    k_min::Float64 = 1.0,
    k_max::Float64 = 10.0,
    r0::Float64 = 0.25,
    eps::Float64 = 0.025,
    m::Int = 2,
    mean::Symbol = :harmonic,
    k_mode::Symbol = :variable,
    k_const::Float64 = 1.0,
    norm::Real = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()

    for nx in nx_list
        ny = nx
        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        r_of_xy = (x, y) -> sqrt((x - center[1])^2 + (y - center[2])^2)
        body = (x, y, _=0) -> r_of_xy(x, y)- R
        capacity = Capacity(body, mesh; method="ImplicitIntegration")
        operator = DiffusionOps(capacity)

        u_exact = (x, y) -> exact_solution(x, y, r_of_xy(x, y), R, m, center)
        f = (x, y, _=0) -> forcing_term(x, y, r_of_xy(x, y), R, m, center; k_min=k_min, k_max=k_max, r0=r0, eps=eps)
        D = if k_mode == :constant
            (x, y, _=0) -> k_const
        else
            (x, y, _=0) -> k_profile(r_of_xy(x, y); k_min=k_min, k_max=k_max, r0=r0, eps=eps)
        end

        bc_outer = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left   => bc_outer,
            :right  => bc_outer,
            :top    => bc_outer,
            :bottom => bc_outer
        ))

        phase = Phase(capacity, operator, f, D)
        solver = DiffusionSteadyMonoVariable(phase, bc_b, Dirichlet((x,y,_=0)->u_exact(x,y)); mean=mean)
        solve_DiffusionSteadyMono!(solver; method=Base.:\)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_exact, solver, capacity, norm, relative)

        push!(h_vals, min(lx / nx, ly / ny))
        push!(err_vals, global_err)
        push!(err_full_vals, full_err)
        push!(err_cut_vals, cut_err)
        push!(err_empty_vals, empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        dx = lx / nx
        dy = ly / ny
        coverage_x = ceil(Int, 2 * R / dx)
        coverage_y = ceil(Int, 2 * R / dy)
        push!(inside_cells_by_dim, [coverage_x, coverage_y])
    end

    return (
        h_vals = h_vals,
        err_vals = err_vals,
        err_full_vals = err_full_vals,
        err_cut_vals = err_cut_vals,
        err_empty_vals = err_empty_vals,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        orders = compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals),
        norm = norm
    )
end

function run_parametric_study(
    nx_list::Vector{Int};
    eps_list::Vector{Float64},
    k_max_list::Vector{Float64},
    k_min::Float64 = 1.0,
    mean_list::Vector{Symbol} = [:harmonic, :arithmetic],
    k_mode_list::Vector{Symbol} = [:variable, :constant],
    k_const::Float64 = 1.0,
    center::Tuple{Float64,Float64} = (0.5, 0.5),
    R::Float64 = 0.5,
    r0::Float64 = 0.25,
    m::Int = 2,
    norm::Real = 2
)
    frames = DataFrame[]

    for k_mode in k_mode_list
        for eps in eps_list
            for k_max in k_max_list
                for mean in mean_list
                    if k_mode == :constant
                        eps != eps_list[1] && continue
                        k_max != k_max_list[1] && continue
                        mean != mean_list[1] && continue
                    end
                data = run_disk_convergence(
                    nx_list;
                    lx=1.0,
                    ly=1.0,
                    center=center,
                    R=R,
                    k_min=k_min,
                    k_max=k_max,
                    r0=r0,
                    eps=eps,
                    m=m,
                    mean=mean,
                    k_mode=k_mode,
                    k_const=k_const,
                    norm=norm,
                    relative=false
                )
                df = make_convergence_dataframe("DiskSteepK", data)
                df.k_min .= k_min
                df.k_max .= k_mode == :constant ? k_const : k_max
                df.k_ratio .= df.k_max ./ k_min
                df.eps .= k_mode == :constant ? missing : eps
                df.r0 .= r0
                df.mean .= string(mean)
                df.k_mode .= string(k_mode)
                df.m .= m
                push!(frames, df) 
                end
            end
        end
    end

    return vcat(frames...)
end

function main(; csv_path=nothing)
    nx_list = [40, 80, 160, 320]
    eps_list = [0.05, 0.025, 0.01]
    k_max_list = [1.0, 10.0, 100.0, 1000.0]

    df = run_parametric_study(
        nx_list;
        eps_list=eps_list,
        k_max_list=k_max_list,
        k_min=1.0,
        mean_list=[:harmonic, :arithmetic],
        k_mode_list=[:variable, :constant],
        k_const=1.0,
        center=(0.5, 0.5),
        R=0.5,
        r0=0.25,
        m=2,
        norm=2
    )

    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "scalar") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ?
        joinpath(results_dir, "Scalar_2D_Diffusion_Poisson_DiskSteepK_Parametric.csv") :
        csv_path
    CSV.write(csv_out, df)

    return (csv_path = csv_out, table = df)
end

results = main()

@testset "Disk steep-k parametric study" begin
    @test isfile(results.csv_path)
    @test nrow(results.table) > 0
end
