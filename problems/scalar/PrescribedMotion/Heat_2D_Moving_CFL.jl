using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using SpecialFunctions
using Roots
using CSV
using Test

"""
Scalar 2D heat equation with an oscillating circular interface. The radius
oscillates periodically, and the exact solution is imposed on the moving
boundary. This script performs mesh convergence studies using CFL-based
interface-velocity time stepping and writes one CSV summary per CFL value.
# Might need to adjust interface centroid computation for moving bodies : bary_interface vs compute_interface_centroid()
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

function oscillating_radius(t, r_mean, r_amp, period)
    return r_mean + r_amp * sin(2π * t / period)
end

function oscillating_velocity(t, r_amp, period)
    return r_amp * (2π / period) * cos(2π * t / period)
end

max_interface_speed(r_amp, period) = abs(r_amp) * (2π / period)

function oscillating_body(center, r_mean, r_amp, period)
    return (x,y,t)-> sqrt((x-center[1])^2 + (y-center[2])^2) - oscillating_radius(t, r_mean, r_amp, period)
end

function run_moving_heat_convergence(
    nx_list::Vector{Int},
    ny_list::Vector{Int},
    r_mean::Float64,
    r_amp::Float64,
    period::Float64,
    center::Tuple{Float64,Float64},
    D::Float64,
    Φ_ana::Function,
    source_term::Function,
    cfl::Float64;
    lx::Float64 = 4.0,
    ly::Float64 = 4.0,
    Tend::Float64 = 0.1,
    norm::Real = 2
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()
    interface_velocity = t -> oscillating_velocity(t, r_amp, period)

    for (nx, ny) in zip(nx_list, ny_list)
        mesh = Penguin.Mesh((nx, ny), (lx, ly), (0.0, 0.0))
        Δh_min = min(lx / nx, ly / ny)
        v_max = max_interface_speed(r_amp, period)
        Δt = v_max < 1e-12 ? 0.5 * Δh_min^2 : cfl * Δh_min / v_max
        Tstart = Δt

        body = oscillating_body(center, r_mean, r_amp, period)
        STmesh = Penguin.SpaceTimeMesh(mesh, [0.0, Δt])
        capacity = Capacity(body, STmesh; method="VOFI", integration_method=:vofijul)
        operator = DiffusionOps(capacity)

        bc = Dirichlet(0.0)
        bc_b = BorderConditions(Dict(
            :left => bc,
            :right => bc,
            :top => bc,
            :bottom => bc
        ))

        interface_bc = Dirichlet((x,y,t)->Φ_ana(x,y,t))
        phase = Phase(capacity, operator, source_term, (x,y,z)->D)

        ndofs = (nx + 1)*(ny + 1)
        u0ₒ = similar(fill(0.0, ndofs))
        for j in 1:ny+1
            for i in 1:nx+1
                idx = (j-1)*(nx+1) + i
                u0ₒ[idx] = Φ_ana(mesh.nodes[1][i], mesh.nodes[2][j], Tstart)
            end
        end
        u0ᵧ = zeros(ndofs)
        u0 = vcat(u0ₒ, u0ᵧ)

        solver = MovingDiffusionUnsteadyMono(phase, bc_b, interface_bc, Δt, u0, mesh, "BE")
        solve_MovingDiffusionUnsteadyMono_cfl!(
            solver, phase, body, Δt, Tstart, Tend, bc_b, interface_bc, mesh, "CN",
            interface_velocity, cfl; method=Base.:\, geometry_method="VOFI", integration_method=:vofijul
        )

        R_tend = oscillating_radius(Tend, r_mean, r_amp, period)
        body_tend = (x,y,_=0) -> sqrt((x-center[1])^2 + (y-center[2])^2) - R_tend
        capacity_tend = Capacity(body_tend, mesh; compute_centroids=false)
        Φ_ana_tend = (x,y)->Φ_ana(x,y,Tend)

        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(Φ_ana_tend, solver, capacity_tend, norm)

        push!(h_vals, lx / nx)
        push!(err_vals, global_err)
        push!(err_full_vals, full_err)
        push!(err_cut_vals, cut_err)
        push!(err_empty_vals, empty_err)
        push!(inside_cells, count_inside_cells(capacity_tend))
        Δx = lx / nx
        Δy = ly / ny
        coverage_x = ceil(Int, 2 * R_tend / Δx)
        coverage_y = ceil(Int, 2 * R_tend / Δy)
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

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_convergence_dataframe(method_name, data)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "scalar") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "$(method_name)_Convergence.csv") : csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function cfl_tag(cfl::Real)
    return replace(string(cfl), "." => "p")
end

function main(; csv_dir=nothing, nx_list=nothing, ny_list=nothing, cfl_list=nothing)
    nx_vals = isnothing(nx_list) ? [4, 8, 16, 32, 64, 128] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list
    cfl_vals = isnothing(cfl_list) ? [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125] : cfl_list
    r_mean = 1.0
    r_amp = 1.0
    period = 1.0
    center = (2.0, 2.0)
    D = 1.0
    Tend = 0.1

    Φ_ana = let center=center, r_mean=r_mean, r_amp=r_amp, period=period
        function (x,y,t)
            r = sqrt((x - center[1])^2 + (y - center[2])^2)
            R_t = oscillating_radius(t, r_mean, r_amp, period)
            if r > R_t
                return 0.0
            end
            return (1 + 0.5 * sin(2π * t / period)) * cos(π * x) * cos(π * y)
        end
    end

    source = let center=center, r_mean=r_mean, r_amp=r_amp, period=period, D=D
        function (x,y,z,t)
            r = sqrt((x - center[1])^2 + (y - center[2])^2)
            R_t = oscillating_radius(t, r_mean, r_amp, period)
            if r > R_t
                return 0.0
            end
            term1 = (π / period) * cos(π * x) * cos(π * y) * cos(2π * t / period)
            term2 = 2 * π^2 * D * (1 + 0.5 * sin(2π * t / period)) * cos(π * x) * cos(π * y)
            return term1 + term2
        end
    end

    data_by_cfl = Dict{Float64, Any}()
    csv_paths = Dict{Float64, String}()

    for cfl in cfl_vals
        data = run_moving_heat_convergence(
            nx_vals, ny_vals, r_mean, r_amp, period, center, D, Φ_ana, source, cfl;
            Tend=Tend
        )
        method_name = "Scalar_2D_Diffusion_Heat_Moving_CFL_$(cfl_tag(cfl))"
        csv_path = isnothing(csv_dir) ? nothing : joinpath(csv_dir, "$(method_name)_Convergence.csv")
        csv_info = write_convergence_csv(method_name, data; csv_path=csv_path)
        data_by_cfl[cfl] = data
        csv_paths[cfl] = csv_info.csv_path
    end

    return (data=data_by_cfl, csv_paths=csv_paths, cfl_vals=cfl_vals)
end

results = main()

@testset "Heat 2D moving CFL convergence" begin
    @test length(results.cfl_vals) == length(results.data)
    for cfl in results.cfl_vals
        data = results.data[cfl]
        @test !isnan(data.orders.all)
        @test length(data.h_vals) == length(data.err_vals)
        @test data.h_vals[1] > data.h_vals[end]
        @test minimum(data.err_vals) < maximum(data.err_vals)
        @test isfile(results.csv_paths[cfl])
    end
end
