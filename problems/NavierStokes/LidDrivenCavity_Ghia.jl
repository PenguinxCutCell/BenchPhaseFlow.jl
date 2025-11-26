using Penguin
using LinearAlgebra
using Statistics
using DelimitedFiles
using Printf
using CSV
using DataFrames
using Test

"""
Steady lid-driven cavity benchmark at Re = 1000 using the Navier–Stokes mono solver.

Matches the classical Ghia et al. (1982) configuration: unit square, no-slip on
three walls, lid speed = 1.0. Solves with Picard then Newton refinement, samples
centerline profiles, compares to reference data stored in `ghia/*.ghia`, and
writes CSV outputs to `results/NavierStokes` (no plots).
"""

###########
# Geometry and discretisation
###########
const nx, ny = 128, 128
const Lx, Ly = 1.0, 1.0
const x0, y0 = 0.0, 0.0

mesh_p = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

###########
# Capacities and operators
###########
full_body = (x, y, _=0.0) -> -1.0
capacity_ux = Capacity(full_body, mesh_ux; compute_centroids=false)
capacity_uy = Capacity(full_body, mesh_uy; compute_centroids=false)
capacity_p  = Capacity(full_body, mesh_p;  compute_centroids=false)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

###########
# Boundary conditions
###########
lid_speed = 1.0

ux_top    = Dirichlet((x, y, t=0.0) -> lid_speed)
ux_bottom = Dirichlet((x, y, t=0.0) -> 0.0)
ux_left   = Dirichlet((x, y, t=0.0) -> 0.0)
ux_right  = Dirichlet((x, y, t=0.0) -> 0.0)

uy_zero = Dirichlet((x, y, t=0.0) -> 0.0)

bc_ux = BorderConditions(Dict(
    :left=>ux_left,
    :right=>ux_right,
    :bottom=>ux_bottom,
    :top=>ux_top,
))

bc_uy = BorderConditions(Dict(
    :left=>uy_zero,
    :right=>uy_zero,
    :bottom=>uy_zero,
    :top=>uy_zero,
))

pressure_gauge = PinPressureGauge()
cut_bc = Dirichlet(0.0)

###########
# Fluid properties and forcing (Re = lid_speed * L / ν = 1000)
###########
const ν = 1e-3
const ρ = 1.0
const μ = ν * ρ
fᵤ = (x, y, z=0.0) -> 0.0
fₚ = (x, y, z=0.0) -> 0.0

fluid = Fluid((mesh_ux, mesh_uy),
              (capacity_ux, capacity_uy),
              (operator_ux, operator_uy),
              mesh_p,
              capacity_p,
              operator_p,
              μ, ρ, fᵤ, fₚ)

###########
# Solver setup
###########
nu_x = prod(operator_ux.size)
nu_y = prod(operator_uy.size)
np = prod(operator_p.size)
x0_vec = zeros(2 * (nu_x + nu_y) + np)  # NavierStokesMono state layout

function load_ghia_profile(path)
    data = readdlm(path)
    coords = data[:, 1] .+ 0.5  # stored as (coord - 0.5)
    values = data[:, 2]
    return coords, values
end

nearest_index(vec, val) = clamp(argmin(abs.(vec .- val)), 1, length(vec))

function main(; csv_dir=nothing)
    solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, cut_bc; x0=x0_vec)

    println("=== Lid-driven cavity (Navier–Stokes mono) ===")
    println("Grid: $(nx) × $(ny), ν = $(ν), lid speed = $(lid_speed)")

    _, picard_iters, picard_res = solve_NavierStokesMono_steady!(
        solver;
        tol=1e-6,
        maxiter=200,
        relaxation=0.7,
        nlsolve_method=:picard,
    )
    println("Picard iterations: iters=$(picard_iters), residual=$(picard_res)")

    _, newton_iters, newton_res = solve_NavierStokesMono_steady!(
        solver;
        tol=1e-9,
        maxiter=40,
        nlsolve_method=:newton,
    )
    println("Newton iterations: iters=$(newton_iters), residual=$(newton_res)")

    blocks = Penguin.navierstokes2D_blocks(solver)
    nu_x = blocks.nu_x
    nu_y = blocks.nu_y

    uωx = solver.x[1:nu_x]
    uγx = solver.x[nu_x+1:2nu_x]
    uωy = solver.x[2nu_x+1:2nu_x+nu_y]
    uγy = solver.x[2nu_x+nu_y+1:2*(nu_x+nu_y)]
    pω  = solver.x[2*(nu_x + nu_y)+1:end]

    mass_residual = blocks.div_x_ω * Vector{Float64}(uωx) +
                    blocks.div_x_γ * Vector{Float64}(uγx) +
                    blocks.div_y_ω * Vector{Float64}(uωy) +
                    blocks.div_y_γ * Vector{Float64}(uγy)

    ke = 0.5 * (sum(abs2, uωx) + sum(abs2, uωy))

    xs = mesh_ux.nodes[1]
    ys = mesh_ux.nodes[2]
    Ux = reshape(uωx, (length(xs), length(ys)))
    Uy = reshape(uωy, (length(xs), length(ys)))

    xs_vis = xs[1:end-1]
    ys_vis = ys[1:end-1]

    icol = nearest_index(xs, x0 + Lx / 2)
    irow = nearest_index(ys, y0 + Ly / 2)
    ux_centerline = vec(Ux[icol, 1:length(ys_vis)])
    uy_centerline = vec(Uy[1:length(xs_vis), irow])

    ghia_y, ghia_ux = load_ghia_profile(joinpath(@__DIR__, "ghia", "xprof.ghia"))
    ghia_x, ghia_uy = load_ghia_profile(joinpath(@__DIR__, "ghia", "yprof.ghia"))

    function sample_profile(grid_coords, field_vals, ref_coords)
        sampled = similar(ref_coords)
        for (k, coord) in pairs(ref_coords)
            idx = nearest_index(grid_coords, coord)
            sampled[k] = field_vals[idx]
        end
        return sampled
    end

    ux_sampled = sample_profile(ys_vis, ux_centerline, ghia_y)
    uy_sampled = sample_profile(xs_vis, uy_centerline, ghia_x)

    ux_diff = ux_sampled .- ghia_ux
    uy_diff = uy_sampled .- ghia_uy

    println("Centerline comparisons vs Ghia et al. (1982):")
    println(@sprintf("  vertical line (u_x): Linf=%.3e, L2=%.3e",
                     maximum(abs, ux_diff), sqrt(mean(abs2, ux_diff))))
    println(@sprintf("  horizontal line (u_y): Linf=%.3e, L2=%.3e",
                     maximum(abs, uy_diff), sqrt(mean(abs2, uy_diff))))

    out_dir = isnothing(csv_dir) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : csv_dir
    mkpath(out_dir)

    summary_df = DataFrame(
        picard_iters = picard_iters,
        picard_residual = picard_res,
        newton_iters = newton_iters,
        newton_residual = newton_res,
        div_linf = maximum(abs, mass_residual),
        kinetic_energy = ke,
        max_abs_ux = maximum(abs, uωx),
        max_abs_uy = maximum(abs, uωy),
        min_p = minimum(pω),
        max_p = maximum(pω),
        linf_ux_ghia = maximum(abs, ux_diff),
        l2_ux_ghia = sqrt(mean(abs2, ux_diff)),
        linf_uy_ghia = maximum(abs, uy_diff),
        l2_uy_ghia = sqrt(mean(abs2, uy_diff))
    )
    summary_path = joinpath(out_dir, "LidDrivenCavity_Ghia_Summary.csv")
    CSV.write(summary_path, summary_df)

    ux_prof_df = DataFrame(y = ghia_y, u_num = ux_sampled, u_ref = ghia_ux)
    uy_prof_df = DataFrame(x = ghia_x, v_num = uy_sampled, v_ref = ghia_uy)
    ux_prof_path = joinpath(out_dir, "LidDrivenCavity_Ghia_UxProfile.csv")
    uy_prof_path = joinpath(out_dir, "LidDrivenCavity_Ghia_UyProfile.csv")
    CSV.write(ux_prof_path, ux_prof_df)
    CSV.write(uy_prof_path, uy_prof_df)

    println("Benchmark completed. Results saved to:")
    println("  ", summary_path)
    println("  ", ux_prof_path)
    println("  ", uy_prof_path)

    return (
        summary = summary_df,
        ux_profile = ux_prof_df,
        uy_profile = uy_prof_df,
        summary_path = summary_path,
        ux_prof_path = ux_prof_path,
        uy_prof_path = uy_prof_path
    )
end

results = main()

@testset "Lid-driven cavity Ghia" begin
    @test size(results.summary, 1) == 1
    @test results.summary.div_linf[1] < 1e-2
    @test isfile(results.summary_path)
    @test isfile(results.ux_prof_path)
    @test isfile(results.uy_prof_path)
end
