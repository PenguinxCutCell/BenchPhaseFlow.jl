using Penguin
using LinearAlgebra
using CSV
using DataFrames
using Test
using Printf

"""
2D unsteady Stokes flow with a cylinder translating in a confined channel.

A rigid cylinder of radius ``a`` moves at constant speed ``U`` in the ``+x``
direction between two parallel plates separated by ``2b`` (confinement ratio
``k = a / b``). We prescribe the cylinder motion, integrate the moving Stokes
solver, and extract drag/lift plus the wall-correction factor
``λ(k) = F_x / (μ U)``. Results are written to
`results/NavierStokes/Stokes_TranslatingCylinder.csv`.
"""

###########
# Geometry and motion
###########
const nx = 96
const ny = 32
const b = 1.0                  # half channel height
const radius = 0.1
const channel_length = 6.0
const x0 = -channel_length / 2
const y0 = -b
const center_x0 = -channel_length / 6
const center_y0 = 0.0
const U_translate = 1.0
const Δt = 0.02
const T_end = 0.4
const scheme = :BE
const geometry_method = "VOFI"

println("Reynolds number: Re = $(2 * radius * U_translate / 1.0)")

k = radius / b

body = let cx0=center_x0, cy0=center_y0, U=U_translate
    (x, y, t) -> begin
        cx = cx0 + U * t
        return radius - sqrt((x - cx)^2 + (y - cy0)^2)  # positive inside the cylinder
    end
end

###########
# Meshes
###########
mesh_p = Penguin.Mesh((nx, ny), (channel_length, 2b), (x0, y0))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (channel_length, 2b), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (channel_length, 2b), (x0, y0 - 0.5 * dy))

###########
# Capacities and operators (initialized at t = 0)
###########
capacity_ux = Capacity((x, y, _=0.0) -> body(x, y, 0.0), mesh_ux)
capacity_uy = Capacity((x, y, _=0.0) -> body(x, y, 0.0), mesh_uy)
capacity_p  = Capacity((x, y, _=0.0) -> body(x, y, 0.0), mesh_p)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

###########
# Boundary conditions
###########
u_zero = Dirichlet(0.0)
bc_ux = BorderConditions(Dict(
    :left   => u_zero,
    :right  => u_zero,
    :bottom => u_zero,
    :top    => u_zero
))
bc_uy = BorderConditions(Dict(
    :left   => u_zero,
    :right  => u_zero,
    :bottom => u_zero,
    :top    => u_zero
))
pressure_gauge = PinPressureGauge()

# Prescribed cylinder velocity: U in x, no motion in y
bc_cut = (
    Dirichlet((x, y, t) -> U_translate),
    Dirichlet((x, y, t) -> 0.0),
)

###########
# Fluid properties
###########
μ = 1.0
ρ = 1.0
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
capacity_kwargs = (; method=geometry_method,
                    integration_method=:vofijul,
                    compute_centroids=true)

nu_x = prod(operator_ux.size)
nu_y = prod(operator_uy.size)
np = prod(operator_p.size)
Ntot = 2 * (nu_x + nu_y) + np
x0_vec = zeros(Ntot)

solver = MovingStokesUnsteadyMono(fluid, (bc_ux, bc_uy), pressure_gauge, bc_cut;
                                  scheme=scheme, x0=x0_vec)

println("=== Translating cylinder in a channel (Stokes) ===")
println(@sprintf("Grid: %d × %d, Δt=%.3f, T_end=%.3f, k=%.3f", nx, ny, Δt, T_end, k))
println(@sprintf("Cylinder radius=%.3f, U=%.3f, channel height=%.3f", radius, U_translate, 2b))

###########
# Time integration
###########
times, states = solve_MovingStokesUnsteadyMono!(solver, body, mesh_p,
                                                Δt, 0.0, T_end,
                                                (bc_ux, bc_uy), bc_cut;
                                                scheme=scheme,
                                                method=Base.:\,
                                                geometry_method=geometry_method,
                                                integration_method=capacity_kwargs.integration_method,
                                                compute_centroids=capacity_kwargs.compute_centroids)

println("Completed $(length(times) - 1) time steps")

###########
# Force diagnostics
###########
function final_force_diagnostics(times, solver, body;
                                 capacity_kwargs, scheme)
    length(times) >= 2 || error("Need at least two time samples to form a space-time mesh")
    t_prev, t_last = times[end - 1], times[end]

    STmesh_ux = Penguin.SpaceTimeMesh(mesh_ux, [t_prev, t_last], tag=mesh_p.tag)
    STmesh_uy = Penguin.SpaceTimeMesh(mesh_uy, [t_prev, t_last], tag=mesh_p.tag)
    STmesh_p  = Penguin.SpaceTimeMesh(mesh_p,  [t_prev, t_last], tag=mesh_p.tag)

    capacity_ux_last = Capacity(body, STmesh_ux; capacity_kwargs...)
    capacity_uy_last = Capacity(body, STmesh_uy; capacity_kwargs...)
    capacity_p_last  = Capacity(body, STmesh_p;  capacity_kwargs...)

    operator_ux_last = DiffusionOps(capacity_ux_last)
    operator_uy_last = DiffusionOps(capacity_uy_last)
    operator_p_last  = DiffusionOps(capacity_p_last)

    block_data = Penguin.stokes2D_moving_blocks(solver.fluid,
                                                (operator_ux_last, operator_uy_last),
                                                (capacity_ux_last, capacity_uy_last),
                                                operator_p_last, capacity_p_last,
                                                scheme)

    force_diag = compute_navierstokes_force_diagnostics(solver, block_data)
    body_force = navierstokes_reaction_force_components(force_diag; acting_on=:body)
    pressure_body = .-force_diag.integrated_pressure
    viscous_body = .-force_diag.integrated_viscous
    coeffs = drag_lift_coefficients(force_diag; ρ=ρ,
                                    U_ref=U_translate,
                                    length_ref=2radius,
                                    acting_on=:body)

    return (; force_diag, body_force, pressure_body, viscous_body, coeffs)
end

diagnostics = final_force_diagnostics(times, solver, body;
                                      capacity_kwargs=capacity_kwargs,
                                      scheme=scheme)

λ_numeric = -diagnostics.body_force[1] / (μ * U_translate)  # positive drag factor

println("Forces at t = $(round(times[end]; digits=3))")
println("  Drag  = $(round(diagnostics.body_force[1]; sigdigits=6)) (pressure=$(round(diagnostics.pressure_body[1]; sigdigits=6)), viscous=$(round(diagnostics.viscous_body[1]; sigdigits=6)))")
println("  Lift  = $(round(diagnostics.body_force[2]; sigdigits=6)) (pressure=$(round(diagnostics.pressure_body[2]; sigdigits=6)), viscous=$(round(diagnostics.viscous_body[2]; sigdigits=6)))")
println("  Cd = $(round(diagnostics.coeffs.Cd; sigdigits=6)), Cl = $(round(diagnostics.coeffs.Cl; sigdigits=6))")
println("  λ(k) = F_x / (μ U) = $(round(λ_numeric; sigdigits=6)) with k = $(k)")

###########
# CSV output
###########
function write_results_csv(diagnostics, λ_value; csv_path=nothing)
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "Stokes_TranslatingCylinder.csv") : csv_path

    df = DataFrame(
        k = k,
        radius = radius,
        b = b,
        U = U_translate,
        Fx = diagnostics.body_force[1],
        Fy = diagnostics.body_force[2],
        Cd = diagnostics.coeffs.Cd,
        Cl = diagnostics.coeffs.Cl,
        lambda = λ_value,
        Δt = Δt,
        T_end = T_end,
        nx = nx,
        ny = ny
    )
    CSV.write(csv_file, df)
    return csv_file
end

csv_path = write_results_csv(diagnostics, λ_numeric)

function main()
    return (
        k = k,
        lambda = λ_numeric,
        Cd = diagnostics.coeffs.Cd,
        Cl = diagnostics.coeffs.Cl,
        Fx = diagnostics.body_force[1],
        Fy = diagnostics.body_force[2],
        csv_path = csv_path,
        time_steps = length(times) - 1
    )
end

results = main()

@testset "Translating cylinder Stokes" begin
    @test isfinite(results.lambda)
    @test isfinite(results.Cd)
    @test isfinite(results.Cl)
    @test isfile(results.csv_path)
    @test results.time_steps >= 1
end

println("Saved drag/lift summary to $(csv_path)")
