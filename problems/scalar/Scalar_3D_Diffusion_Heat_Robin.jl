using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using Roots
using CSV
using Test
using CairoMakie

"""
Scalar 3D Diffusion Heat Equation with Robin Boundary Conditions
This benchmark mirrors the 2D Robin problem but inside a sphere. The initial
temperature is uniform (w0) and the interface satisfies ∂ₙw + k w = 0. The
analytical solution uses the eigenvalues μₙ obtained from μ cot(μ) + kR - 1 = 0:
w(r,t) = (2kR² w0)/r * Σ Cₙ sin(μₙ r / R) exp(-a μₙ² t / R²),
where Cₙ = sin(μₙ)[μₙ²+(kR-1)²]/(μₙ²[μₙ²+kR(kR-1)]). Th
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

function robin_mu_roots(N::Int, k::Float64, R::Float64; eps::Float64 = 1e-6)
    eq(mu) = mu * cot(mu) + k * R - 1
    deq(mu) = cot(mu) - mu / sin(mu)^2
    roots = zeros(Float64, N)
    for n in 1:N
        guess = (n - 0.5) * π
        try
            println("Finding root $n with initial guess $guess with Newton-Raphson")
            roots[n] = find_zero((eq, deq), guess, Roots.Newton(); atol=1e-12, rtol=1e-12, maxiters=100)
        catch
            println("  ... Newton-Raphson failed, falling back to bisection")
            left = (n - 1) * π + eps
            right = n * π - eps
            roots[n] = find_zero(eq, (left, right), Roots.Bisection(); atol=1e-12, maxiters=200)
        end
    end
    return roots
end

function robin_spherical_heat_solution(center::Tuple{Float64,Float64,Float64}, R::Float64;
    t::Float64 = 0.1,
    k::Float64 = 1.0,
    a::Float64 = 1.0,
    w0::Float64 = 1.0,
    Nroots::Int = 400,
    tol::Float64 = 1e-12
)
    μs = robin_mu_roots(Nroots, k, R)
    println(μs)
    coeffs = similar(μs)
    for (i, μ) in enumerate(μs)
        coeffs[i] = sin(μ) * (μ^2 + (k*R - 1)^2) / (μ^2 * (μ^2 + k*R*(k*R - 1)))
    end

    prefactor = 2 * k * R^2 * w0

    return function (x::Float64, y::Float64, z::Float64)
        r = sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2)
        if r >= R
            return w0
        end
        if r < 1e-12
            s = 0.0
            for (μ, C) in zip(μs, coeffs)
                term = C * μ / R * exp(-a * μ^2 * t / R^2)
                s += term
                if abs(term) < tol
                    break
                end
            end
            return prefactor * s
        else
            s = 0.0
            for (μ, C) in zip(μs, coeffs)
                term = C * sin(μ * r / R) * exp(-a * μ^2 * t / R^2)
                s += term
                if abs(term) < tol
                    break
                end
            end
            return prefactor * s / r
        end
    end
end

function run_heat_robin_convergence_3d(
    nx_list::Vector{Int},
    ny_list::Vector{Int},
    nz_list::Vector{Int},
    radius::Float64,
    center::Tuple{Float64,Float64,Float64},
    u_analytical::Function;
    lx::Float64 = 4.0,
    ly::Float64 = 4.0,
    lz::Float64 = 4.0,
    norm::Int = 2,
    Tend::Float64 = 0.1,
    k::Float64 = 1.0,
    w0::Float64 = 1.0
)
    h_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()
    last_solver = nothing
    last_capacity = nothing

    for (nx, ny, nz) in zip(nx_list, ny_list, nz_list)
        mesh = Penguin.Mesh((nx, ny, nz), (lx, ly, lz), (0.0, 0.0, 0.0))
        sphere = (x,y,z,_=0) -> sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2) - radius
        capacity = Capacity(sphere, mesh)
        operator = DiffusionOps(capacity)

        bc_boundary = Robin(k, 1.0, 0.0)
        bc_b = BorderConditions(Dict(
            :left   => Dirichlet(1.0),
            :right  => Dirichlet(1.0),
            :top    => Dirichlet(1.0),
            :bottom => Dirichlet(1.0),
            :front  => Dirichlet(1.0),
            :back   => Dirichlet(1.0)
        ))
        phase = Phase(capacity, operator, (x,y,z,t)->0.0, (x,y,z)->1.0)

        ndofs = (nx + 1) * (ny + 1) * (nz + 1)
        u0ₒ = ones(ndofs) * w0
        u0ᵧ = zeros(ndofs)
        u0 = vcat(u0ₒ, u0ᵧ)

        Δx = lx / nx
        Δy = ly / ny
        Δz = lz / nz
        Δt = 0.5 * min(Δx^2, Δy^2, Δz^2)

        solver = DiffusionUnsteadyMono(phase, bc_b, bc_boundary, Δt, u0, "BE")
        solve_DiffusionUnsteadyMono!(solver, phase, Δt, Tend, bc_b, bc_boundary, "CN"; method=Base.:\)
        _, _, global_err, full_err, cut_err, empty_err =
            check_convergence(u_analytical, solver, capacity, norm)

        push!(h_vals, min(Δx, Δy, Δz))
        push!(err_vals, global_err)
        push!(err_full_vals, full_err)
        push!(err_cut_vals, cut_err)
        push!(err_empty_vals, empty_err)
        push!(inside_cells, count_inside_cells(capacity))
        coverage_x = ceil(Int, 2 * radius / Δx)
        coverage_y = ceil(Int, 2 * radius / Δy)
        coverage_z = ceil(Int, 2 * radius / Δz)
        push!(inside_cells_by_dim, [coverage_x, coverage_y, coverage_z])
        last_solver = solver
        last_capacity = capacity
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
        norm = norm,
        last_solver = last_solver,
        last_capacity = last_capacity
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

function plot_radial_profile_heat_robin(u_analytical::Function, solver, capacity;
    center::Tuple{Float64,Float64,Float64} = (2.0, 2.0, 2.0),
    radius::Float64 = 1.0,
    output_path::Union{Nothing,String} = nothing,
    n_reference::Int = 400
)
    u_ana, u_num, _, _, _, _ = check_convergence(u_analytical, solver, capacity, 2)
    mesh = capacity.mesh
    nodes = mesh.nodes
    dy = nodes[2][2] - nodes[2][1]
    dz = nodes[3][2] - nodes[3][1]
    line_tol = 0.5 * min(dy, dz)
    radial_tol = radius + line_tol

    r_vals = Float64[]
    num_vals = Float64[]
    ana_vals = Float64[]
    line_indices = Int[]
    for (idx, (c, num, ana)) in enumerate(zip(capacity.C_ω, u_num, u_ana))
        r = sqrt((c[1] - center[1])^2 + (c[2] - center[2])^2 + (c[3] - center[3])^2)
        if r <= radial_tol && abs(c[2] - center[2]) <= line_tol && abs(c[3] - center[3]) <= line_tol && c[1] >= center[1]
            push!(r_vals, r)
            push!(num_vals, num)
            push!(ana_vals, ana)
            push!(line_indices, idx)
        end
    end

    isempty(r_vals) && error("No grid points found on the radial line y=$(center[2]), z=$(center[3]).")
    perm = sortperm(r_vals)
    r_line = r_vals[perm]
    num_line = num_vals[perm]
    ana_line = ana_vals[perm]
    line_indices = line_indices[perm]

    r_plot = range(0, stop=radius, length=n_reference)
    ana_plot = [u_analytical(center[1] + r, center[2], center[3]) for r in r_plot]
    # all values >radius are replace by NaN in the analytical solution
    ana_plot .= ifelse.(r_plot .< radius, ana_plot, NaN)

    u_num_gamma = solver.x[end÷2+1:end]
    length(u_num_gamma) == length(u_num) || @warn "Unexpected gamma length ($(length(u_num_gamma))) vs bulk ($(length(u_num))); using overlap for plotting."
    valid_mask = line_indices .<= length(u_num_gamma)
    gamma_r = r_line[valid_mask]
    gamma_vals = u_num_gamma[line_indices[valid_mask]]
    if !isempty(gamma_r)
        idx_max = argmax(gamma_r)
        gamma_vals = [gamma_vals[idx_max]]
        gamma_r = [radius]
    end

    fig = Figure(resolution=(720, 420))
    ax = Axis(fig[1, 1],
        title = "Scalar 3D Robin : radial profile",
        xlabel = "r (distance from center)",
        ylabel = "Φ value",)
    lines!(ax, r_plot, ana_plot, label="Analytical", color=:black)
    scatter!(ax, r_line, num_line, label="Numerical (y=z=$(center[2]))", color=:dodgerblue, markersize=7)
    if !isempty(gamma_r)
        scatter!(ax, gamma_r, gamma_vals, label="Interface value", color=:red, markersize=8)
    end
    axislegend(ax, position=:rt)
    xlims!(ax, 0, radius+0.1)

    profile_path = isnothing(output_path) ? joinpath(BENCH_ROOT, "results", "scalar", "Scalar_3D_Diffusion_Heat_Robin_Profile.png") : output_path
    mkpath(dirname(profile_path))
    save(profile_path, fig)
    return (
        path = profile_path,
        r = r_line,
        numerical = num_line,
        analytical = ana_line,
        gamma_r = gamma_r,
        gamma_vals = gamma_vals,
        figure = fig
    )
end

function main(; csv_path=nothing, nx_list=nothing, ny_list=nothing, nz_list=nothing, plot_profile::Bool=false, profile_path=nothing)
    nx_vals = isnothing(nx_list) ? [2, 4, 8, 16, 32, 64, 96] : nx_list
    ny_vals = isnothing(ny_list) ? nx_vals : ny_list
    nz_vals = isnothing(nz_list) ? nx_vals : nz_list
    radius = 1.0
    center = (2.0, 2.0, 2.0)
    Tend = 0.1
    k = 1.0
    u_analytical = robin_spherical_heat_solution(center, radius; t=Tend, k=k, a=1.0, w0=1.0, Nroots=200)

    data = run_heat_robin_convergence_3d(
        nx_vals, ny_vals, nz_vals, radius, center, u_analytical;
        lx = 4.0, ly = 4.0, lz = 4.0, norm = 2, Tend = Tend
    )

    csv_info = write_convergence_csv("Scalar_3D_Diffusion_Heat_Robin", data; csv_path=csv_path)
    profile = nothing
    if plot_profile
        data.last_solver === nothing && error("No solver available for plotting; ensure at least one mesh run.")
        data.last_capacity === nothing && error("No capacity available for plotting; ensure at least one mesh run.")
        profile = plot_radial_profile_heat_robin(
            u_analytical, data.last_solver, data.last_capacity;
            center = center, radius = radius, output_path = profile_path
        )
    end
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table, profile = profile)
end

results = main(plot_profile=true)

@testset "Heat 3D Robin convergence" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test orders.all > 1.0
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
    println("orders: ", orders)
end
