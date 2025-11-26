using Penguin
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using Statistics
using CSV
using DataFrames
using Test

"""
2D Navier–Stokes Poiseuille benchmark (steady) with cut-cell channel (no plots).

Adapted from `examples/2D/Stokes/poiseuille_2d_cut.jl` but uses the steady
Navier–Stokes mono solver. Computes mid-column profile errors against the
analytical parabola and writes results to `results/NavierStokes`.
"""

const nx, ny = 64, 64
const Lx, Ly = 2.0, 1.0
const x0, y0 = 0.0, 0.0

const y_wall_bot = 0.2
const y_wall_top = 0.8
const channel_height = y_wall_top - y_wall_bot

const Umax = 1.0
const μ = 1.0
const ρ = 1.0

body = (x, y, _=0) -> begin
    if y < y_wall_bot
        return y_wall_bot - y
    elseif y > y_wall_top
        return y - y_wall_top
    else
        return -min(y - y_wall_bot, y_wall_top - y)
    end
end

mesh_p = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
dx, dy = Lx / nx, Ly / ny
mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

capacity_ux = Capacity(body, mesh_ux; compute_centroids=false)
capacity_uy = Capacity(body, mesh_uy; compute_centroids=false)
capacity_p  = Capacity(body, mesh_p; compute_centroids=false)
operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

parabola = (x, y) -> begin
    if y < y_wall_bot || y > y_wall_top
        return 0.0
    else
        y_local = y - y_wall_bot
        return 4 * Umax * y_local * (channel_height - y_local) / (channel_height^2)
    end
end

ux_left  = Dirichlet(parabola)
ux_right = Dirichlet(parabola)
ux_bot   = Dirichlet((x, y)->0.0)
ux_top   = Dirichlet((x, y)->0.0)
bc_ux = BorderConditions(Dict(
    :left=>ux_left, :right=>ux_right, :bottom=>ux_bot, :top=>ux_top
))

uy_zero = Dirichlet((x, y)->0.0)
bc_uy = BorderConditions(Dict(:left=>uy_zero, :right=>uy_zero, :bottom=>uy_zero, :top=>uy_zero))

pressure_gauge = PinPressureGauge()
u_bc = Dirichlet(0.0)

fᵤ = (x, y, z=0.0) -> 0.0
fₚ = (x, y, z=0.0) -> 0.0

fluid = Fluid((mesh_ux, mesh_uy),
              (capacity_ux, capacity_uy),
              (operator_ux, operator_uy),
              mesh_p,
              capacity_p,
              operator_p,
              μ, ρ, fᵤ, fₚ)

nu = prod(operator_ux.size)
np = prod(operator_p.size)
x0_vec = zeros(4 * nu + np)

function run_poiseuille_cut()
    solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, u_bc; x0=x0_vec)
    _, iters, res = solve_NavierStokesMono_steady!(solver; tol=1e-10, maxiter=200, relaxation=1.0)
    println("Navier–Stokes steady solve finished: iterations=", iters, ", residual=", res)

    uωx = solver.x[1:nu]
    uωy = solver.x[2nu+1:3nu]
    pω  = solver.x[4nu+1:end]

    xs = mesh_ux.nodes[1]
    ys = mesh_ux.nodes[2] .+ 0.5 * dy  # mid-column y-locations for ux
    LIux = LinearIndices((length(xs), length(ys)))

    icol = Int(cld(length(xs), 2))
    ux_profile = [uωx[LIux[icol, j]] for j in 1:length(ys)]
    ux_analytical = [parabola(0.0, y) for y in ys]

    fluid_indices = findall(y -> y_wall_bot <= y <= y_wall_top, ys)
    profile_err = ux_profile[fluid_indices] .- ux_analytical[fluid_indices]
    ℓ2_profile = sqrt(sum(abs2, profile_err) / length(profile_err))
    ℓinf_profile = maximum(abs, profile_err[2:end-1])  # ignore boundary points

    j_center = argmin(abs.(ys .- (y_wall_bot + y_wall_top) / 2))
    ux_centerline = ux_profile[j_center]
    ux_analytical_centerline = parabola(0.0, ys[j_center])
    rel_center_err = abs(ux_centerline - ux_analytical_centerline) /
        (abs(ux_analytical_centerline) + eps())

    println("Profile L2 error (fluid region) = ", ℓ2_profile, ", Linf = ", ℓinf_profile)
    println("Relative centerline error = ", rel_center_err)
    println("checks:")
    println("- profile Linf < 1e-2? ", ℓinf_profile < 1e-2)
    println("- centerline relative error < 1e-2? ", rel_center_err < 1e-2)
    println("- residual < 1e-10? ", res < 1e-10)

    return (
        iterations = iters,
        residual = res,
        l2_profile = ℓ2_profile,
        linf_profile = ℓinf_profile,
        rel_center_err = rel_center_err,
        mean_ux = mean(ux_profile),
        max_ux = maximum(abs.(ux_profile)),
        mean_uy = mean(uωy),
        max_uy = maximum(abs.(uωy)),
        mean_p = mean(pω),
        min_p = minimum(pω),
        max_p = maximum(pω)
    )
end

function write_results_csv(results; csv_path=nothing)
    df = DataFrame([results])
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "Poiseuille_2D_Cut.csv") : csv_path
    CSV.write(csv_file, df)
    return csv_file
end

function main(; csv_path=nothing)
    results = run_poiseuille_cut()
    csv_file = write_results_csv(results; csv_path=csv_path)
    return (results=results, csv_path=csv_file)
end

results = main()

@testset "Poiseuille 2D cut Navier–Stokes" begin
    @test isfinite(results.results.l2_profile)
    @test results.results.residual < 1e-6
    @test isfile(results.csv_path)
end

println("Poiseuille 2D cut benchmark completed. Results saved to ", results.csv_path)
