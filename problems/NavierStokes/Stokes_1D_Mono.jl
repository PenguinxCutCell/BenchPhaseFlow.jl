using Penguin
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using Statistics
using CSV
using DataFrames
using Test

"""
1D monophasic Stokes benchmark (no plotting).

This is a CSV-only version of `examples/1D/Stokes/stokes_mono.jl`. It solves a
staggered 1D Stokes system with homogeneous velocity Dirichlet boundaries and a
pinned pressure gauge, then writes key diagnostics to
`results/NavierStokes/Stokes_1D_Mono.csv`.
"""

struct Stokes1DParams
    nx::Int
    lx::Float64
    x0::Float64
    xint::Float64
    μ::Float64
    ρ::Float64
end

Stokes1DParams(; nx=160, lx=4.0, x0=0.0, xint=2.0, μ=1.0, ρ=1.0) =
    Stokes1DParams(nx, lx, x0, xint, μ, ρ)

function run_stokes_1d(params::Stokes1DParams)
    nx, lx, x0, xint, μ, ρ = params.nx, params.lx, params.x0, params.xint, params.μ, params.ρ

    mesh_p = Penguin.Mesh((nx,), (lx,), (x0,))
    dx = lx / nx
    mesh_u = Penguin.Mesh((nx,), (lx,), (x0 - 0.5 * dx,))

    body = (x, _=0) -> (x - xint)^2 - 1.0
    capacity_u = Capacity(body, mesh_u)
    capacity_p = Capacity(body, mesh_p)
    operator_u = DiffusionOps(capacity_u)
    operator_p = DiffusionOps(capacity_p)

    u_left = Dirichlet(0.0)
    u_right = Dirichlet(0.0)
    bc_u = BorderConditions(Dict(:bottom => u_left, :top => u_right))
    pressure_gauge = PinPressureGauge()
    u_bc = Dirichlet(1.0)

    fᵤ = (x, y=0.0, z=0.0) -> 0.0
    fₚ = (x, y=0.0, z=0.0) -> 0.0

    fluid = Fluid(mesh_u, capacity_u, operator_u, mesh_p, capacity_p, operator_p, μ, ρ, fᵤ, fₚ)

    nu = prod(operator_u.size)
    np = prod(operator_p.size)
    u0ₒ = zeros(nu)
    u0ᵧ = zeros(nu)
    p0ₒ = zeros(np)
    x0_vec = vcat(u0ₒ, u0ᵧ, p0ₒ)

    solver = StokesMono(fluid, bc_u, pressure_gauge, u_bc; x0=x0_vec)
    solve_StokesMono!(solver)

    uω = solver.x[1:nu]
    uγ = solver.x[nu+1:2nu]
    pω = solver.x[2nu+1:end]

    mask_u = diag(capacity_u.V) .> 1e-12
    mask_p = diag(capacity_p.V) .> 1e-12
    uω = uω[mask_u]
    uγ = uγ[mask_u]
    pω = pω[mask_p]

    mean_uω = mean(uω)
    max_uω = maximum(abs.(uω))
    mean_uγ = mean(uγ)
    max_uγ = maximum(abs.(uγ))
    min_p = minimum(pω)
    max_p = maximum(pω)
    mean_p = mean(pω)

    uω_rand = randn(nu)
    uγ_zero = zeros(nu)
    pω_rand = randn(np)
    div_u = - (operator_p.G' + operator_p.H') * uω_rand + (operator_p.H') * uγ_zero
    grad_p = operator_p.Wꜝ * (operator_p.G + operator_p.H) * pω_rand
    lhs = dot(operator_p.V * div_u, pω_rand)
    rhs = dot(uω_rand, grad_p)
    adj_rel = abs(lhs - rhs) / max(1.0, abs(lhs), abs(rhs))

    Iμ⁻¹ = build_I_D(operator_u, 1 / μ, capacity_u)
    WG = operator_u.Wꜝ * operator_u.G
    S = Iμ⁻¹ * (operator_u.G' * WG)
    sym_err = opnorm(Matrix(S - S')) / max(1e-12, opnorm(Matrix(S)))
    uω_t = randn(nu)
    rq = dot(uω_t, S * uω_t)

    return (
        mean_uω = mean_uω,
        max_uω = max_uω,
        mean_uγ = mean_uγ,
        max_uγ = max_uγ,
        min_p = min_p,
        max_p = max_p,
        mean_p = mean_p,
        adjoint_rel_error = adj_rel,
        viscous_sym_error = sym_err,
        rayleigh_quotient = rq
    )
end

function write_results_csv(results; csv_path=nothing)
    df = DataFrame([results]) # wrap NamedTuple as a single row table
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "Stokes_1D_Mono.csv") : csv_path
    CSV.write(csv_file, df)
    return csv_file
end

function main(; params::Stokes1DParams=Stokes1DParams(), csv_path=nothing)
    results = run_stokes_1d(params)
    csv_file = write_results_csv(results; csv_path=csv_path)
    return (results=results, csv_path=csv_file)
end

results = main()

@testset "Stokes 1D mono benchmark" begin
    res = results.results
    @test isfinite(res.mean_uω)
    @test isfinite(res.mean_p)
    @test res.rayleigh_quotient >= -1e-10
    @test isfile(results.csv_path)
end

println("Stokes 1D mono benchmark completed. Results saved to ", results.csv_path)
