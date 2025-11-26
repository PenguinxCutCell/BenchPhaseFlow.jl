using Penguin
using Statistics
using Printf
using CSV
using DataFrames
using Test

"""
Steady Navier–Stokes flow around a circular cylinder: wake length benchmark (CSV only).

Solves steady Navier–Stokes for subcritical Re values, measures wake length
behind the cylinder (reattachment point on centerline), and compares to simple
literature brackets. Outputs summary CSV to `results/NavierStokes`.
"""

###########
# Geometry and discretisation
###########
const nx, ny = 256, 128
const channel_length = 4.0
const channel_height = 1.0
const x0, y0 = -0.5, -0.5

const circle_center = (0.5, 0.0)
const circle_radius = 0.2
const diameter = 2circle_radius

circle_body = (x, y, _=0.0) -> circle_radius - hypot(x - circle_center[1], y - circle_center[2])

mesh_p  = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0, y0))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0, y0 - 0.5 * dy))

capacity_ux = Capacity(circle_body, mesh_ux)
capacity_uy = Capacity(circle_body, mesh_uy)
capacity_p  = Capacity(circle_body, mesh_p)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

###########
# Boundary conditions
###########
const U∞ = 1.0
uniform_inlet = (x, y, t=0.0) -> U∞

ux_left   = Dirichlet((x, y, t=0.0) -> uniform_inlet(x, y, t))
ux_right  = Outflow()
ux_bottom = Symmetry()
ux_top    = Symmetry()

uy_zero   = Dirichlet((x, y, t=0.0) -> 0.0)

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
interface_bc = Dirichlet(0.0)

###########
# Solver helpers
###########
const ρ = 1.0
fᵤ = (x, y, z=0.0, t=0.0) -> 0.0
fₚ = (x, y, z=0.0, t=0.0) -> 0.0

nu_x = prod(operator_ux.size)
nu_y = prod(operator_uy.size)
np = prod(operator_p.size)
x_template = zeros(2 * (nu_x + nu_y) + np)

function build_fluid(μ::Float64)
    Fluid((mesh_ux, mesh_uy),
          (capacity_ux, capacity_uy),
          (operator_ux, operator_uy),
          mesh_p,
          capacity_p,
          operator_p,
          μ, ρ, fᵤ, fₚ)
end

function solve_steady_state(μ; picard_tol=1e-9, picard_maxiter=120, relaxation=1.0)
    fluid = build_fluid(μ)
    solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, interface_bc; x0=copy(x_template))

    _, picard_iters, picard_res = solve_NavierStokesMono_steady!(
        solver; nlsolve_method=:picard, tol=picard_tol, maxiter=picard_maxiter, relaxation=relaxation)

    blocks = Penguin.navierstokes2D_blocks(solver)
    mass_residual = blocks.div_x_ω * Vector{Float64}(solver.x[1:nu_x]) +
                    blocks.div_x_γ * Vector{Float64}(solver.x[nu_x+1:2nu_x]) +
                    blocks.div_y_ω * Vector{Float64}(solver.x[2nu_x+1:2nu_x+nu_y]) +
                    blocks.div_y_γ * Vector{Float64}(solver.x[2nu_x+nu_y+1:2*(nu_x+nu_y)])

    return solver, picard_iters, picard_res, maximum(abs, mass_residual)
end

###########
# Wake-length utilities
###########
xs = mesh_ux.nodes[1]
ys = mesh_ux.nodes[2]
x_trailing = circle_center[1] + circle_radius

nearest_index(vec, val) = clamp(argmin(abs.(vec .- val)), 1, length(vec))

function recirculation_length(Ux::Matrix{Float64})
    j_center = nearest_index(ys, circle_center[2])
    start_idx = searchsortedfirst(xs, x_trailing + dx)
    start_idx > length(xs) && return (0.0, missing)

    first_neg = nothing
    for i in start_idx:length(xs)
        if Ux[i, j_center] < 0
            first_neg = i
            break
        end
    end
    first_neg === nothing && return (0.0, missing)

    last_neg = first_neg
    next_pos = nothing
    for i in (first_neg + 1):length(xs)
        if Ux[i, j_center] < 0
            last_neg = i
        else
            next_pos = i
            break
        end
    end

    next_pos === nothing && return (NaN, missing)

    x_neg = xs[last_neg]
    x_pos = xs[next_pos]
    u_neg = Ux[last_neg, j_center]
    u_pos = Ux[next_pos, j_center]
    x_zero = x_neg - u_neg * (x_pos - x_neg) / (u_pos - u_neg)

    wake_len = max(0.0, x_zero - x_trailing)
    return (wake_len, x_zero)
end

###########
# Reference data (simple brackets)
###########
wake_brackets = Dict(
    20 => (0.73, 0.94),
    40 => (1.60, 2.29),
)

###########
# Benchmark sweep
###########
function run_wake_benchmark(Re_list::Vector{Int})
    results = Vector{NamedTuple}(undef, length(Re_list))
    println("=== Steady flow around a circular cylinder benchmark ===")

    for (i, Re) in enumerate(Re_list)
        μ = ρ * U∞ * diameter / Re
        println()
        println(@sprintf("Re = %d  -> viscosity μ = %.4e", Re, μ))

        solver, picard_iters, picard_res, mass_inf = solve_steady_state(μ)
        println(@sprintf("  Picard:  iters=%d  residual=%.3e", picard_iters, picard_res))
        println(@sprintf("  max|div u| ≈ %.3e", mass_inf))

        Ux = reshape(solver.x[1:nu_x], (length(xs), length(ys)))
        uy_offset = 2nu_x
        Uy = reshape(solver.x[uy_offset+1:uy_offset+nu_y], (length(xs), length(ys)))

        wake_len, x_reattach = recirculation_length(Ux)
        wake_len_D = wake_len / diameter
        println(@sprintf("  wake length L_r = %.3e (%.3f D)", wake_len, wake_len_D))

        bracket = get(wake_brackets, Re, (NaN, NaN))
        within_bracket = isfinite(bracket[1]) ? (bracket[1] <= wake_len_D <= bracket[2]) : false
        if isfinite(bracket[1])
            println(@sprintf("  bracket check: [%.2f, %.2f] D -> %s",
                             bracket[1], bracket[2], within_bracket ? "inside" : "outside"))
        else
            println("  no wake bracket available")
        end

        results[i] = (
            Re = Re,
            viscosity = μ,
            wake_length = wake_len,
            wake_length_D = wake_len_D,
            reattach_x = x_reattach,
            mass_inf = mass_inf,
            picard_iters = picard_iters,
            bracket_min = bracket[1],
            bracket_max = bracket[2],
            within_bracket = within_bracket
        )
    end

    return results
end

function write_results_csv(results; csv_path=nothing)
    df = DataFrame(results)
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "CylinderWake_Steady.csv") : csv_path
    CSV.write(csv_file, df)
    return (csv_file=csv_file, table=df)
end

function main(; Re_list=[20, 40], csv_path=nothing)
    res = run_wake_benchmark(Re_list)
    csv_info = write_results_csv(res; csv_path=csv_path)
    return (results=res, csv_path=csv_info.csv_file, table=csv_info.table)
end

results = main()

@testset "Cylinder wake steady benchmark" begin
    @test isfile(results.csv_path)
    @test length(results.results) == length([20, 40])
    @test all(isfinite.(results.table.wake_length_D))
end

println("Cylinder wake steady benchmark complete. Results saved to ", results.csv_path)
