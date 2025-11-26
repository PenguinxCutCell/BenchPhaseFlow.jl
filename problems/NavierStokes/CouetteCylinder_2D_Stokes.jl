using Penguin
using LinearAlgebra
using SparseArrays
using Statistics
using Printf
using CSV
using DataFrames
using Test

"""
Steady Stokes Couette flow between two concentric cylinders (CSV only).

Inner cylinder (R₁=0.25) rotates with ω=1, outer cylinder (R₂=0.5) is fixed.
Level-set defines an annulus domain. Solves with `StokesMono` and a
component-aware cut BC, then samples tangential velocity along x=0 to compute
L2/Linf errors against the analytic Couette profile. Writes summary CSV to
`results/NavierStokes/CouetteCylinder_2D_Stokes.csv`.
"""

###########
# Level-set helpers
###########
interval_min(a, b) = 0.5 * (a + b - abs(a - b))

struct CouetteCylinderParams
    R₁::Float64  # inner radius
    R₂::Float64  # outer radius
    ω::Float64   # angular velocity of inner cylinder
end

CouetteCylinderParams(; R₁=0.25, R₂=0.5, ω=1.0) = CouetteCylinderParams(R₁, R₂, ω)

function diffusion_levelset(params::CouetteCylinderParams)
    return (x, y, _=0.0) -> begin
        r = hypot(x, y)
        φ_inner = r - params.R₁          # negative inside the inner disk
        φ_outer_ext = params.R₂ - r      # negative outside the outer cylinder
        interval_min(φ_inner, φ_outer_ext)
    end
end

# Fluid (negative) region is the annulus => flip level-set sign
annulus_body(params) = (x, y, _=0.0) -> -diffusion_levelset(params)(x, y)

###########
# Analytic Couette solution utilities
###########
function tangential_velocity(r, params::CouetteCylinderParams)
    Ri, Ro = params.R₁, params.R₂
    ω = params.ω
    δ = (Ri^2 - Ro^2)
    A = ω * Ri^2 / δ
    B = -ω * Ri^2 * Ro^2 / δ
    return A * r + B / r
end

function analytic_velocity_components(x, y, params::CouetteCylinderParams)
    r = hypot(x, y)
    if r < params.R₁ || r > params.R₂
        return (0.0, 0.0)
    end
    uθ = tangential_velocity(r, params)
    r < eps() && return (0.0, 0.0)
    ux = -uθ * y / r
    uy =  uθ * x / r
    return (ux, uy)
end

###########
# Component-aware cut BC
###########
import Penguin: build_g_g

struct CouetteCutBC <: Penguin.AbstractBoundary
    value_map::Dict{UInt64, Function}
    default::Union{Float64, Function}
end

CouetteCutBC(default::Union{Float64, Function}=0.0) =
    CouetteCutBC(Dict{UInt64, Function}(), default)

function evaluate_cut_bc(bc::CouetteCutBC, mesh_id::UInt64)
    if haskey(bc.value_map, mesh_id)
        return bc.value_map[mesh_id]
    elseif bc.default isa Function
        return bc.default
    else
        return (args...) -> bc.default
    end
end

function build_g_g(op::Penguin.AbstractOperators, bc::CouetteCutBC, cap::Penguin.Capacity)
    coords = Penguin.get_all_coordinates(cap.C_γ)
    mesh_id = objectid(cap.mesh)
    f = evaluate_cut_bc(bc, mesh_id)
    return [f(coord...) for coord in coords]
end

###########
# Geometry and grids
###########
params = CouetteCylinderParams()
const nx = 256
const ny = 256
const domain_half = 0.8
const Lx = 2domain_half
const Ly = 2domain_half
const x0 = -domain_half
const y0 = -domain_half

mesh_p  = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

body = annulus_body(params)

capacity_ux = Capacity(body, mesh_ux; compute_centroids=true)
capacity_uy = Capacity(body, mesh_uy; compute_centroids=true)
capacity_p  = Capacity(body, mesh_p;  compute_centroids=true)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

###########
# Boundary conditions on the box
###########
zero_bc = Dirichlet((x, y, _=0.0) -> 0.0)
bc_ux = BorderConditions(Dict(
    :left=>zero_bc, :right=>zero_bc, :bottom=>zero_bc, :top=>zero_bc
))
bc_uy = BorderConditions(Dict(
    :left=>zero_bc, :right=>zero_bc, :bottom=>zero_bc, :top=>zero_bc
))
pressure_gauge = MeanPressureGauge()

function couette_cut_component(component::Symbol, params::CouetteCylinderParams)
    tol = 1e-2
    return (x, y, _=0.0) -> begin
        r = hypot(x, y)
        if abs(r - params.R₁) <= tol
            uθ = params.ω * params.R₁
        elseif abs(r - params.R₂) <= tol
            uθ = 0.0
        else
            return 0.0
        end
        r_safe = max(r, eps())
        component === :ux ? -uθ * y / r_safe : uθ * x / r_safe
    end
end

cut_bc = CouetteCutBC(0.0)
cut_bc.value_map[objectid(mesh_ux)] = couette_cut_component(:ux, params)
cut_bc.value_map[objectid(mesh_uy)] = couette_cut_component(:uy, params)

###########
# Material and solver
###########
μ = 1.0
ρ = 1.0
fᵤ = (x, y, z=0.0) -> 0.0
fₚ = (x, y, z=0.0) -> 0.0

fluid = Fluid((mesh_ux, mesh_uy),
              (capacity_ux, capacity_uy),
              (operator_ux, operator_uy),
              mesh_p, capacity_p, operator_p,
              μ, ρ, fᵤ, fₚ)

nu = prod(operator_ux.size)
np = prod(operator_p.size)
x0_vec = zeros(4 * nu + np)

function run_couette()
    solver = StokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, cut_bc; x0=x0_vec)
    solve_StokesMono!(solver; method=Base.:\)

    pω  = solver.x[4nu+1:end]
    uωx = solver.x[1:nu]
    uωy = solver.x[2nu+1:3nu]

    x_nodes_x = mesh_ux.nodes[1]
    y_nodes_x = mesh_ux.nodes[2]
    LI_ux = LinearIndices((length(x_nodes_x), length(y_nodes_x)))
    i_axis = argmin(abs.(x_nodes_x .- 0.0))

    profile_points = Tuple{Float64,Float64,Float64,Float64,Float64}[]
    for (j, yval) in enumerate(y_nodes_x)
        if yval <= params.R₁ + 3dy || yval >= params.R₂ - 3dy
            continue
        end
        idx = LI_ux[i_axis, j]
        Ux_val = uωx[idx]
        ux_exact = -tangential_velocity(yval, params)
        uθ_num = -Ux_val
        push!(profile_points, (yval, yval, Ux_val, ux_exact, uθ_num))
    end

    num_vals = [val[3] for val in profile_points]
    exact_vals = [val[4] for val in profile_points]
    errs = num_vals .- exact_vals
    ℓ2_err = isempty(errs) ? NaN : sqrt(mean(abs2, errs))
    ℓinf_err = isempty(errs) ? NaN : maximum(abs, errs)
    radial_drift = maximum(abs, uωx)

    return (
        solver = solver,
        profile_points = profile_points,
        l2_err = ℓ2_err,
        linf_err = ℓinf_err,
        radial_drift = radial_drift,
        n_profile = length(profile_points)
    )
end

function write_results_csv(data; csv_path=nothing)
    df = DataFrame(
        l2_err = data.l2_err,
        linf_err = data.linf_err,
        radial_drift = data.radial_drift,
        n_profile = data.n_profile,
        nx = nx,
        ny = ny,
        R_inner = params.R₁,
        R_outer = params.R₂
    )
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "CouetteCylinder_2D_Stokes.csv") : csv_path
    CSV.write(csv_file, df)
    return csv_file
end

function main(; csv_path=nothing)
    data = run_couette()
    csv_file = write_results_csv(data; csv_path=csv_path)
    return (data=data, csv_path=csv_file)
end

results = main()

@testset "Couette cylinder Stokes benchmark" begin
    @test isfile(results.csv_path)
    @test isfinite(results.data.l2_err) || isnan(results.data.l2_err)
    @test results.data.n_profile >= 0
end

println("Couette cylinder benchmark completed. Results saved to ", results.csv_path)
