using Penguin
using LinearAlgebra
using Printf
using CSV
using DataFrames
using Test

"""
Taylor–Green vortex convergence for the Navier–Stokes prototype (no plotting).

Solves the manufactured Taylor–Green flow on [0, 2π]² with exact velocity
Dirichlet boundaries, computes weighted L2 errors for u and v, estimates
pairwise convergence rates, and writes a CSV to `results/NavierStokes`.
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
const scheme = :CN  # Crank–Nicolson for viscous terms

u_exact(x, y, t) =  sin(k * x) * cos(k * y) * exp(-2.0 * ν * k^2 * t)
v_exact(x, y, t) = -cos(k * x) * sin(k * y) * exp(-2.0 * ν * k^2 * t)

interior_indices(n) = 2:(n - 1)  # skip boundaries for error norms

"""
Run the convergence study for a list of resolutions `ns`.
Returns a NamedTuple of vectors (h, err_u, err_v, rate_u, rate_v).
"""
function run_taylor_green(ns::Vector{Int})
    errors_u = Float64[]
    errors_v = Float64[]
    hs = Float64[]

    for n in ns
        nx = n; ny = n
        mesh_p  = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
        dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
        dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
        Δt = 0.25 * min(dx, dy)
        mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
        mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

        body = (x, y, _=0.0) -> -1.0  # whole domain is fluid

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
        bc_cut = Dirichlet(0.0)

        fᵤ = (x, y, z=0.0) -> 0.0
        fₚ = (x, y, z=0.0) -> 0.0

        fluid = Fluid((mesh_ux, mesh_uy),
                      (cap_ux, cap_uy),
                      (op_ux, op_uy),
                      mesh_p,
                      cap_p,
                      op_p,
                      μ, ρ, fᵤ, fₚ)

        nu_x = prod(op_ux.size)
        nu_y = prod(op_uy.size)
        np = prod(op_p.size)
        Ntot = 2 * (nu_x + nu_y) + np
        x0_vec = zeros(Float64, Ntot)

        xs_ux = mesh_ux.nodes[1]; ys_ux = mesh_ux.nodes[2]
        xs_uy = mesh_uy.nodes[1]; ys_uy = mesh_uy.nodes[2]

        u_init = zeros(Float64, nu_x)
        v_init = zeros(Float64, nu_y)

        idx = 1
        for j in 1:length(ys_ux), i in 1:length(xs_ux)
            u_init[idx] = u_exact(xs_ux[i], ys_ux[j], 0.0)
            idx += 1
        end

        idx = 1
        for j in 1:length(ys_uy), i in 1:length(xs_uy)
            v_init[idx] = v_exact(xs_uy[i], ys_uy[j], 0.0)
            idx += 1
        end

        x0_vec[1:nu_x] .= u_init
        x0_vec[nu_x+1:2nu_x] .= u_init          # tie DOFs
        x0_vec[2nu_x+1:2nu_x+nu_y] .= v_init
        x0_vec[2nu_x+nu_y+1:2*(nu_x+nu_y)] .= v_init
        x0_vec[2*(nu_x+nu_y)+1:end] .= 0.0

        solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, bc_cut; x0=x0_vec)

        @printf("Running Navier–Stokes resolution %d×%d ...\n", nx, ny)
        _, states = solve_NavierStokesMono_unsteady!(solver; Δt=Δt, T_end=t_end, scheme=scheme, method=Base.:\)
        final = states[end]

        u_num = final[1:nu_x]
        v_num = final[2nu_x+1:2nu_x+nu_y]

        Ux = reshape(u_num, length(xs_ux), length(ys_ux))
        Uy = reshape(v_num, length(xs_uy), length(ys_uy))

        Ux_ex = [u_exact(xs_ux[i], ys_ux[j], t_end) for i in 1:length(xs_ux), j in 1:length(ys_ux)]
        Uy_ex = [v_exact(xs_uy[i], ys_uy[j], t_end) for i in 1:length(xs_uy), j in 1:length(ys_uy)]

        ix_range_u = interior_indices(length(xs_ux))
        iy_range_u = interior_indices(length(ys_ux))
        ix_range_v = interior_indices(length(xs_uy))
        iy_range_v = interior_indices(length(ys_uy))

        Vux = diag(op_ux.V)
        Vuy = diag(op_uy.V)

        function weighted_L2_grid(num, ex, mask_i, mask_j, Vdiag; remove_mean=false)
            ni, nj = size(num, 1), size(num, 2)
            total_w = 0.0
            weighted_sum = 0.0
            for j in 1:nj, i in 1:ni
                if (i in mask_i) && (j in mask_j)
                    lin = (j - 1) * ni + i
                    w = Vdiag[lin]
                    err = num[i, j] - ex[i, j]
                    total_w += w
                    weighted_sum += w * err
                end
            end
            mean_err = (remove_mean && total_w > 0) ? weighted_sum / total_w : 0.0

            accum = 0.0
            for j in 1:nj, i in 1:ni
                if (i in mask_i) && (j in mask_j)
                    lin = (j - 1) * ni + i
                    w = Vdiag[lin]
                    err = (num[i, j] - ex[i, j]) - mean_err
                    accum += w * err^2
                end
            end
            return sqrt(accum)
        end

        err_u = weighted_L2_grid(Ux, reshape(Ux_ex, size(Ux)), ix_range_u, iy_range_u, Vux)
        err_v = weighted_L2_grid(Uy, reshape(Uy_ex, size(Uy)), ix_range_v, iy_range_v, Vuy)

        push!(errors_u, err_u)
        push!(errors_v, err_v)
        push!(hs, max(Lx / nx, Ly / ny))

        @printf("  h=%.5e  ||u||_L2=%.5e  ||v||_L2=%.5e\n", hs[end], err_u, err_v)
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

    return (h=hs, err_u=errors_u, err_v=errors_v, rate_u=r_u, rate_v=r_v)
end

function write_results_csv(data; csv_path=nothing)
    n = length(data.h)
    rates_u = vcat(NaN, data.rate_u...)
    rates_v = vcat(NaN, data.rate_v...)
    df = DataFrame(
        h = data.h,
        error_u = data.err_u,
        error_v = data.err_v,
        rate_u = rates_u,
        rate_v = rates_v
    )
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "TaylorGreen_2D_Mono.csv") : csv_path
    CSV.write(csv_file, df)
    return csv_file
end

function main(; ns=[8, 16, 32, 64, 128, 256], csv_path=nothing)
    data = run_taylor_green(ns)
    csv_file = write_results_csv(data; csv_path=csv_path)
    return (data=data, csv_path=csv_file)
end

results = main()

@testset "Taylor–Green 2D mono convergence" begin
    @test length(results.data.h) == length(results.data.err_u) == length(results.data.err_v)
    @test results.data.h[1] > results.data.h[end]
    @test minimum(results.data.err_u) < maximum(results.data.err_u)
    @test minimum(results.data.err_v) < maximum(results.data.err_v)
    @test isfile(results.csv_path)
end

println("Taylor–Green convergence completed. Results saved to ", results.csv_path)
