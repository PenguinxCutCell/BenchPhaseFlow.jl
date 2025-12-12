using Penguin
using LinearAlgebra
using Printf
using CSV
using DataFrames
using Test

"""
Taylor–Green vortex with an embedded circular body (no plotting).

Solves the manufactured Taylor–Green flow on [0, 2π]² with a circular cut-cell
obstacle centered in the domain. Exact Dirichlet data is imposed on the domain
borders and on the embedded boundary. Computes weighted L2 errors for u and v,
estimates convergence rates, and writes a CSV to `results/NavierStokes`.
"""

const Lx = 2π
const Ly = 2π
const x0 = 0.0
const y0 = 0.0

const μ = 1.0
const ρ = 1.0
const k = 1.0
const ν = μ / ρ

const t_end = 0.1
const Δt = 0.001
const scheme = :CN  # Crank–Nicolson for viscous terms

const circle_center = (Lx / 2, Ly / 2)
const circle_radius = Lx / 6

u_exact(x, y, t) =  sin(k * x) * cos(k * y) * exp(-2.0 * ν * k^2 * t)
v_exact(x, y, t) = -cos(k * x) * sin(k * y) * exp(-2.0 * ν * k^2 * t)

"""
Evaluate an analytical field on capacity centroids for any supported dimension.
"""
function eval_on_centroids(capacity::Capacity{D}, f::Function) where D
    centroids = capacity.C_ω
    if D == 1
        return map(c -> f(c[1]), centroids)
    elseif D == 2
        return map(c -> f(c[1], c[2]), centroids)
    elseif D == 3
        return map(c -> f(c[1], c[2], c[3]), centroids)
    elseif D == 4
        return map(c -> f(c[1], c[2], c[3], c[4]), centroids)
    else
        error("Unsupported dimension: $D")
    end
end

"""
Boolean mask for points that sit on the outer domain boundary (used to drop
boundary values from error norms).
"""
function domain_border_mask(capacity::Capacity{2})
    tol = 1e-12
    mask = falses(length(capacity.C_ω))
    x_min, x_max = x0, x0 + Lx
    y_min, y_max = y0, y0 + Ly
    for (i, c) in pairs(capacity.C_ω)
        x, y = c
        if abs(x - x_min) ≤ tol || abs(x - x_max) ≤ tol ||
           abs(y - y_min) ≤ tol || abs(y - y_max) ≤ tol
            mask[i] = true
        end
    end
    return mask
end

function component_norms(err, capacity::Capacity, u_ana, p::Real, relative::Bool)
    cell_types = capacity.cell_types
    border_mask = domain_border_mask(capacity)
    drop_border = i -> !border_mask[i]

    idx_all = filter(drop_border, findall((cell_types .== 1) .| (cell_types .== -1)))
    idx_full = filter(drop_border, findall(cell_types .== 1))
    idx_cut = filter(drop_border, findall(cell_types .== -1))
    idx_empty = filter(drop_border, findall(cell_types .== 0))

    if relative
        global_err = Penguin.relative_lp_norm(err, idx_all, p, capacity, u_ana)
        full_err = Penguin.relative_lp_norm(err, idx_full, p, capacity, u_ana)
        cut_err = Penguin.relative_lp_norm(err, idx_cut, p, capacity, u_ana)
        empty_err = Penguin.relative_lp_norm(err, idx_empty, p, capacity, u_ana)
    else
        global_err = Penguin.lp_norm(err, idx_all, p, capacity)
        full_err = Penguin.lp_norm(err, idx_full, p, capacity)
        cut_err = Penguin.lp_norm(err, idx_cut, p, capacity)
        empty_err = Penguin.lp_norm(err, idx_empty, p, capacity)
    end

    return (globall=global_err, full=full_err, cut=cut_err, empty=empty_err)
end

"""
    velocity_convergence_metrics(u_analytical::Function, v_analytical::Function,
                                 solver, cap_ux::Capacity{D}, cap_uy::Capacity{D};
                                 p::Real=2, relative::Bool=false) where D

Compute weighted Lᵖ errors for the x- and y-velocity components stored in
`solver.x` with layout `[uωx, uγx, uωy, uγy, pω]`.
"""
function velocity_convergence_metrics(u_analytical::Function, v_analytical::Function,
                                      solver, cap_ux::Capacity{D}, cap_uy::Capacity{D};
                                      p::Real=2, relative::Bool=false) where D
    u_ana = eval_on_centroids(cap_ux, u_analytical)
    v_ana = eval_on_centroids(cap_uy, v_analytical)

    nu_x = length(u_ana)
    nu_y = length(v_ana)

    uωx = solver.x[1:nu_x]
    uωy = solver.x[2nu_x + 1:2nu_x + nu_y]

    err_x = u_ana .- uωx
    err_y = v_ana .- uωy

    ux_norms = component_norms(err_x, cap_ux, u_ana, p, relative)
    uy_norms = component_norms(err_y, cap_uy, v_ana, p, relative)

    println("u_x L$p norms -> all=$(ux_norms.globall), full=$(ux_norms.full), cut=$(ux_norms.cut), empty=$(ux_norms.empty)")
    println("u_y L$p norms -> all=$(uy_norms.globall), full=$(uy_norms.full), cut=$(uy_norms.cut), empty=$(uy_norms.empty)")

    return (
        u_ana=u_ana,
        v_ana=v_ana,
        u_num=uωx,
        v_num=uωy,
        ux_errors=ux_norms,
        uy_errors=uy_norms
    )
end

"""
Run the convergence study for a list of resolutions `ns`.
Returns a NamedTuple of vectors with mesh size, errors, and rates.
"""
function run_taylor_green_circle(ns::Vector{Int})
    errors_u = Float64[]
    errors_v = Float64[]
    cut_errors_u = Float64[]
    cut_errors_v = Float64[]
    hs = Float64[]

    for n in ns
        nx = n; ny = n
        mesh_p  = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
        dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
        dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
        mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
        mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

        body = (x, y, _=0.0) -> -1.0 #-(circle_radius - sqrt((x - circle_center[1])^2 + (y - circle_center[2])^2))

        cap_ux = Capacity(body, mesh_ux)
        cap_uy = Capacity(body, mesh_uy)
        cap_p  = Capacity(body, mesh_p)

        op_ux = DiffusionOps(cap_ux)
        op_uy = DiffusionOps(cap_uy)
        op_p  = DiffusionOps(cap_p)

        bc_ux = BorderConditions(Dict(
            :left   => Dirichlet((x, y, t)->u_exact(x, y, t)),
            :right  => Dirichlet((x, y, t)->u_exact(x, y, t)),
            :bottom => Dirichlet((x, y, t)->u_exact(x, y, t)),
            :top    => Dirichlet((x, y, t)->u_exact(x, y, t)),
        ))

        bc_uy = BorderConditions(Dict(
            :left   => Dirichlet((x, y, t)->v_exact(x, y, t)),
            :right  => Dirichlet((x, y, t)->v_exact(x, y, t)),
            :bottom => Dirichlet((x, y, t)->v_exact(x, y, t)),
            :top    => Dirichlet((x, y, t)->v_exact(x, y, t)),
        ))

        pressure_gauge = PinPressureGauge()
        bc_cut = Dirichlet((x, y, t)->u_exact(x, y, t))  # applied to both interface velocity components

        fᵤ = (x, y, z=0.0) -> 0.0
        fₚ = (x, y, z=0.0) -> 0.0

        fluid = Fluid((mesh_ux, mesh_uy),
                      (cap_ux, cap_uy),
                      (op_ux, op_uy),
                      mesh_p,
                      cap_p,
                      op_p,
                      μ, ρ, fᵤ, fₚ)

        nu = prod(op_ux.size)
        @assert nu == prod(op_uy.size) "Taylor–Green circle expects matching staggered DOF counts"
        np = prod(op_p.size)
        Ntot = 4 * nu + np
        x0_vec = zeros(Float64, Ntot)

        xs_ux = mesh_ux.nodes[1]; ys_ux = mesh_ux.nodes[2]
        xs_uy = mesh_uy.nodes[1]; ys_uy = mesh_uy.nodes[2]

        u_init = [u_exact(x, y, 0.0) for x in xs_ux, y in ys_ux] |> vec
        v_init = [v_exact(x, y, 0.0) for x in xs_uy, y in ys_uy] |> vec

        x0_vec[1:nu] .= u_init
        x0_vec[nu+1:2nu] .= u_init
        x0_vec[2nu+1:3nu] .= v_init
        x0_vec[3nu+1:4nu] .= v_init
        x0_vec[4nu+1:end] .= 0.0

        solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, bc_cut; x0=x0_vec)

        @printf("Running Navier–Stokes Taylor–Green with circle, resolution %d×%d ...\n", nx, ny)
        _, states = solve_NavierStokesMono_unsteady!(solver; Δt=Δt, T_end=t_end, scheme=scheme, method=Base.:\)
        solver.x .= states[end]

        metrics = velocity_convergence_metrics(
            (x, y) -> u_exact(x, y, t_end),
            (x, y) -> v_exact(x, y, t_end),
            solver, cap_ux, cap_uy;
            p=2, relative=false
        )

        push!(errors_u, metrics.ux_errors.global)
        push!(errors_v, metrics.uy_errors.global)
        push!(cut_errors_u, metrics.ux_errors.cut)
        push!(cut_errors_v, metrics.uy_errors.cut)
        push!(hs, max(Lx / nx, Ly / ny))

        @printf("  h=%.5e  ||u||_L2=%.5e  ||v||_L2=%.5e  (cut u=%.5e, cut v=%.5e)\n",
                hs[end], errors_u[end], errors_v[end], cut_errors_u[end], cut_errors_v[end])
    end

    function rate(h, e)
        r = Float64[]
        for i in 2:length(e)
            push!(r, log(e[i] / e[i-1]) / log(h[i] / h[i-1]))
        end
        return r
    end

    r_u = rate(hs, errors_u)
    r_v = rate(hs, errors_v)

    return (h=hs, err_u=errors_u, err_v=errors_v,
            err_u_cut=cut_errors_u, err_v_cut=cut_errors_v,
            rate_u=r_u, rate_v=r_v)
end

function write_results_csv(data; csv_path=nothing)
    n = length(data.h)
    rates_u = vcat(NaN, data.rate_u...)
    rates_v = vcat(NaN, data.rate_v...)
    df = DataFrame(
        h = data.h,
        error_u = data.err_u,
        error_v = data.err_v,
        error_u_cut = data.err_u_cut,
        error_v_cut = data.err_v_cut,
        rate_u = rates_u,
        rate_v = rates_v
    )
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "TaylorGreen_2D_Circle.csv") : csv_path
    CSV.write(csv_file, df)
    return csv_file
end

function main(; ns=[4, 8, 16, 32, 64], csv_path=nothing)
    data = run_taylor_green_circle(ns)
    csv_file = write_results_csv(data; csv_path=csv_path)
    return (data=data, csv_path=csv_file)
end

results = main()

@testset "Taylor–Green 2D circle convergence" begin
    @test length(results.data.h) == length(results.data.err_u) == length(results.data.err_v)
    @test results.data.h[1] > results.data.h[end]
    @test minimum(results.data.err_u) < maximum(results.data.err_u)
    @test minimum(results.data.err_v) < maximum(results.data.err_v)
    @test isfile(results.csv_path)
end

println("Taylor–Green convergence with circular body completed. Results saved to ", results.csv_path)
