using Penguin
using LinearAlgebra
using Printf
using CSV
using DataFrames
using Test

"""
3D steady Stokes flow around a spherical obstacle (no plots).

The flow is driven by a uniform free-stream velocity in the x-direction and a
no-slip sphere embedded via a level-set capacity. Boundary velocities are
imposed from the analytical creeping-flow solution around a sphere so the
numerical solution should converge to the exact field. Results are written to
`results/NavierStokes/Stokes_3D_Sphere.csv`.
"""

const Lx = 10.0
const Ly = 10.0
const Lz = 10.0
const x0 = -Lx / 2
const y0 = -Ly / 2
const z0 = -Lz / 2

const sphere_center = (0.0, 0.0, 0.0)
const sphere_radius = 1.0
const U_inf = 1.0
const μ = 1.0
const ρ = 1.0

sphere_body = (x, y, z, _=0.0) -> sphere_radius - sqrt((x - sphere_center[1])^2 +
                                                       (y - sphere_center[2])^2 +
                                                       (z - sphere_center[3])^2)

function stokes_sphere_velocity(x, y, z; U=U_inf, a=sphere_radius, center=sphere_center)
    dx = x - center[1]
    dy = y - center[2]
    dz = z - center[3]
    r = sqrt(dx^2 + dy^2 + dz^2)

    if r <= a
        return (0.0, 0.0, 0.0)
    end

    cos_theta = dx / r
    sin_theta = sqrt(dy^2 + dz^2) / r

    vr = U * (1.0 - 3.0 * a / (2.0 * r) + a^3 / (2.0 * r^3)) * cos_theta
    vtheta = -U * (1.0 - 3.0 * a / (4.0 * r) - a^3 / (4.0 * r^3)) * sin_theta

    if sin_theta > 1e-12
        rho = sqrt(dy^2 + dz^2)
        vx = vr * cos_theta - vtheta * sin_theta
        vy = vr * dy / r + vtheta * cos_theta * dy / rho
        vz = vr * dz / r + vtheta * cos_theta * dz / rho
    else
        vx = vr * cos_theta
        vy = 0.0
        vz = 0.0
    end

    return (vx, vy, vz)
end

function weighted_L2_interior(num::AbstractVector, exact::AbstractVector, weights::AbstractVector, dims::NTuple{3,Int})
    nx, ny, nz = dims
    num3 = reshape(num, dims)
    ex3 = reshape(exact, dims)
    w3 = reshape(weights, dims)

    ix = nx > 4 ? (2:(nx - 2)) : (1:nx)
    iy = ny > 4 ? (2:(ny - 2)) : (1:ny)
    iz = nz > 4 ? (2:(nz - 2)) : (1:nz)

    acc = 0.0
    wsum = 0.0
    @inbounds for k in iz, j in iy, i in ix
        w = w3[i, j, k]
        err = num3[i, j, k] - ex3[i, j, k]
        acc += w * err^2
        wsum += w
    end
    return wsum > 0 ? sqrt(acc / wsum) : NaN
end

function pairwise_rates(hs::Vector{Float64}, errs::Vector{Float64})
    n = length(hs)
    rates = fill(NaN, n)
    for i in 2:n
        if hs[i-1] != hs[i] && errs[i-1] > 0 && errs[i] > 0
            rates[i] = log(errs[i-1] / errs[i]) / log(hs[i-1] / hs[i])
        end
    end
    return rates
end

function run_stokes_sphere(ns::Vector{Int})
    results = NamedTuple[]

    println("=== 3D Stokes flow around a sphere: convergence ===")
    println(@sprintf("%6s %8s %12s %12s %12s", "N", "h", "||u_x||2", "||u_y||2", "||u_z||2"))

    for N in ns
        mesh_p = Penguin.Mesh((N, N, N), (Lx, Ly, Lz), (x0, y0, z0))
        Δx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
        Δy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
        Δz = mesh_p.nodes[3][2] - mesh_p.nodes[3][1]

        mesh_ux = Penguin.Mesh((N, N, N), (Lx, Ly, Lz), (x0 - 0.5 * Δx, y0, z0))
        mesh_uy = Penguin.Mesh((N, N, N), (Lx, Ly, Lz), (x0, y0 - 0.5 * Δy, z0))
        mesh_uz = Penguin.Mesh((N, N, N), (Lx, Ly, Lz), (x0, y0, z0 - 0.5 * Δz))

        cap_ux = Capacity(sphere_body, mesh_ux; method="VOFI", integration_method=:vofijul)
        cap_uy = Capacity(sphere_body, mesh_uy; method="VOFI", integration_method=:vofijul)
        cap_uz = Capacity(sphere_body, mesh_uz; method="VOFI", integration_method=:vofijul)
        cap_p = Capacity(sphere_body, mesh_p; method="VOFI", integration_method=:vofijul)

        op_ux = DiffusionOps(cap_ux)
        op_uy = DiffusionOps(cap_uy)
        op_uz = DiffusionOps(cap_uz)
        op_p = DiffusionOps(cap_p)

        bc_ux = BorderConditions(Dict(
            :left   => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[1]),
            :right  => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[1]),
            :bottom => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[1]),
            :top    => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[1]),
            :front  => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[1]),
            :back   => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[1]),
        ))
        bc_uy = BorderConditions(Dict(
            :left   => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[2]),
            :right  => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[2]),
            :bottom => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[2]),
            :top    => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[2]),
            :front  => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[2]),
            :back   => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[2]),
        ))
        bc_uz = BorderConditions(Dict(
            :left   => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[3]),
            :right  => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[3]),
            :bottom => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[3]),
            :top    => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[3]),
            :front  => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[3]),
            :back   => Dirichlet((x, y, z, t=0.0) -> stokes_sphere_velocity(x, y, z)[3]),
        ))

        fᵤ = (x, y, z, t=0.0) -> 0.0
        fₚ = (x, y, z, t=0.0) -> 0.0

        fluid = Fluid((mesh_ux, mesh_uy, mesh_uz),
                      (cap_ux, cap_uy, cap_uz),
                      (op_ux, op_uy, op_uz),
                      mesh_p,
                      cap_p,
                      op_p,
                      μ, ρ, fᵤ, fₚ)

        nu_x = prod(op_ux.size)
        nu_y = prod(op_uy.size)
        nu_z = prod(op_uz.size)
        np = prod(op_p.size)
        total_dofs = 2 * (nu_x + nu_y + nu_z) + np

        xs_ux, ys_ux, zs_ux = mesh_ux.nodes
        xs_uy, ys_uy, zs_uy = mesh_uy.nodes
        xs_uz, ys_uz, zs_uz = mesh_uz.nodes

        Ux_exact = [stokes_sphere_velocity(x, y, z)[1] for x in xs_ux, y in ys_ux, z in zs_ux]
        Uy_exact = [stokes_sphere_velocity(x, y, z)[2] for x in xs_uy, y in ys_uy, z in zs_uy]
        Uz_exact = [stokes_sphere_velocity(x, y, z)[3] for x in xs_uz, y in ys_uz, z in zs_uz]

        x0_vec = zeros(Float64, total_dofs)
        x0_vec[1:nu_x] .= vec(Ux_exact)
        x0_vec[nu_x + 1:2*nu_x] .= vec(Ux_exact)
        x0_vec[2*nu_x + 1:2*nu_x + nu_y] .= vec(Uy_exact)
        x0_vec[2*nu_x + nu_y + 1:2*nu_x + 2*nu_y] .= vec(Uy_exact)
        x0_vec[2*nu_x + 2*nu_y + 1:2*nu_x + 2*nu_y + nu_z] .= vec(Uz_exact)
        x0_vec[2*nu_x + 2*nu_y + nu_z + 1:2 * (nu_x + nu_y + nu_z)] .= vec(Uz_exact)

        pressure_gauge = PinPressureGauge()
        bc_cut = Dirichlet(0.0)

        solver = StokesMono(fluid, (bc_ux, bc_uy, bc_uz), pressure_gauge, bc_cut; x0=x0_vec)
        solve_StokesMono!(solver)

        uωx = solver.x[1:nu_x]
        uωy = solver.x[2*nu_x + 1:2*nu_x + nu_y]
        uωz = solver.x[2*nu_x + 2*nu_y + 1:2*nu_x + 2*nu_y + nu_z]

        err_ux = weighted_L2_interior(uωx, vec(Ux_exact), diag(op_ux.V), Tuple(op_ux.size))
        err_uy = weighted_L2_interior(uωy, vec(Uy_exact), diag(op_uy.V), Tuple(op_uy.size))
        err_uz = weighted_L2_interior(uωz, vec(Uz_exact), diag(op_uz.V), Tuple(op_uz.size))

        h = maximum((Lx / N, Ly / N, Lz / N))
        push!(results, (; N, h, err_ux, err_uy, err_uz))

        println(@sprintf("%6d %8.4f %12.4e %12.4e %12.4e", N, h, err_ux, err_uy, err_uz))
    end

    return results
end

function write_results_csv(results; csv_path=nothing)
    hs = [r.h for r in results]
    err_ux = [r.err_ux for r in results]
    err_uy = [r.err_uy for r in results]
    err_uz = [r.err_uz for r in results]

    rate_ux = pairwise_rates(hs, err_ux)
    rate_uy = pairwise_rates(hs, err_uy)
    rate_uz = pairwise_rates(hs, err_uz)

    df = DataFrame(
        N = [r.N for r in results],
        h = hs,
        err_ux = err_ux,
        err_uy = err_uy,
        err_uz = err_uz,
        rate_ux = rate_ux,
        rate_uy = rate_uy,
        rate_uz = rate_uz
    )

    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "Stokes_3D_Sphere.csv") : csv_path
    CSV.write(csv_file, df)
    return (csv_path=csv_file, rates=(ux=rate_ux, uy=rate_uy, uz=rate_uz), table=df)
end

function main(; ns=[8, 12, 16, 24, 32], csv_path=nothing)
    res = run_stokes_sphere(ns)
    csv_info = write_results_csv(res; csv_path=csv_path)
    return (results=res, csv_path=csv_info.csv_path, table=csv_info.table, rates=csv_info.rates)
end

results = main()

@testset "Stokes 3D sphere convergence" begin
    @test length(results.results) == length(results.table.h)
    @test results.results[1].h > results.results[end].h
    @test minimum([r.err_ux for r in results.results]) < maximum([r.err_ux for r in results.results])
    @test minimum([r.err_uy for r in results.results]) < maximum([r.err_uy for r in results.results])
    @test minimum([r.err_uz for r in results.results]) < maximum([r.err_uz for r in results.results])
    @test isfile(results.csv_path)
    @test count(!isnan, results.rates.ux) >= 1
end

println("Stokes 3D sphere convergence completed. Results saved to ", results.csv_path)
