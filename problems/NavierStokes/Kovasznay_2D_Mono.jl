using Penguin
using Statistics
using Printf
using LinearAlgebra
using CSV
using DataFrames
using Test

"""
Kovasznay flow convergence study (steady Navier–Stokes, Re ≈ 40) — CSV only.

Solves the manufactured Kovasznay solution on [0,1]² for a set of uniform grids,
enforcing exact Dirichlet data on all boundaries. Reports L2/Linf errors for
u, v, p and estimated pairwise convergence rates. Outputs are written to
`results/NavierStokes/Kovasznay_2D_Mono.csv` (no plots).
"""

const Re = 40.0
const ρ = 1.0
const ν = 1.0 / Re
const μ = ρ * ν
const λ = Re / 2 - sqrt(Re^2 / 4 + 4π^2)

kovasznay_u(x, y) = 1 - exp(λ * x) * cos(2π * y)
kovasznay_v(x, y) = (λ / (2π)) * exp(λ * x) * sin(2π * y)
kovasznay_p(x, y) = 0.5 * (1 - exp(2λ * x))

trim_copy(A) = size(A, 1) > 1 && size(A, 2) > 1 ? copy(@view A[1:end-1, 1:end-1]) : copy(A)
trim_coords(v) = length(v) > 1 ? copy(v[1:end-1]) : copy(v)

function run_kovasznay(resolutions::Vector{Int})
    results = NamedTuple[]

    println("=== Kovasznay flow steady convergence (Re ≈ 40) ===")
    println(@sprintf("%6s %7s %12s %12s %12s %12s %12s %12s",
                    "N", "h", "‖u‖₂ err", "‖u‖∞ err", "‖v‖₂ err", "‖v‖∞ err",
                    "‖p‖₂ err", "‖p‖∞ err"))

    for N in resolutions
        nx = N
        ny = N
        Lx, Ly = 1.0, 1.0
        x0, y0 = 0.0, 0.0

        mesh_p  = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0))
        dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
        dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
        mesh_ux = Penguin.Mesh((nx, ny), (Lx, Ly), (x0 - 0.5 * dx, y0))
        mesh_uy = Penguin.Mesh((nx, ny), (Lx, Ly), (x0, y0 - 0.5 * dy))

        body = (x, y, _=0.0) -> -1.0
        capacity_ux = Capacity(body, mesh_ux)
        capacity_uy = Capacity(body, mesh_uy)
        capacity_p  = Capacity(body, mesh_p)

        operator_ux = DiffusionOps(capacity_ux)
        operator_uy = DiffusionOps(capacity_uy)
        operator_p  = DiffusionOps(capacity_p)

        bc_ux = BorderConditions(Dict(
            :left   => Dirichlet((x, y, t=0.0) -> kovasznay_u(x, y)),
            :right  => Dirichlet((x, y, t=0.0) -> kovasznay_u(x, y)),
            :bottom => Dirichlet((x, y, t=0.0) -> kovasznay_u(x, y)),
            :top    => Dirichlet((x, y, t=0.0) -> kovasznay_u(x, y)),
        ))
        bc_uy = BorderConditions(Dict(
            :left   => Dirichlet((x, y, t=0.0) -> kovasznay_v(x, y)),
            :right  => Dirichlet((x, y, t=0.0) -> kovasznay_v(x, y)),
            :bottom => Dirichlet((x, y, t=0.0) -> kovasznay_v(x, y)),
            :top    => Dirichlet((x, y, t=0.0) -> kovasznay_v(x, y)),
        ))

        fᵤ = (x, y, z=0.0, t=0.0) -> 0.0
        fₚ = (x, y, z=0.0, t=0.0) -> 0.0

        fluid = Fluid((mesh_ux, mesh_uy),
                      (capacity_ux, capacity_uy),
                      (operator_ux, operator_uy),
                      mesh_p,
                      capacity_p,
                      operator_p,
                      μ, ρ, fᵤ, fₚ)

        nu_x = prod(operator_ux.size)
        nu_y = prod(operator_uy.size)
        np = prod(operator_p.size)
        total_dofs = 2 * (nu_x + nu_y) + np

        xs_ux = mesh_ux.nodes[1]
        ys_ux = mesh_ux.nodes[2]
        xs_uy = mesh_uy.nodes[1]
        ys_uy = mesh_uy.nodes[2]
        xs_p = mesh_p.nodes[1]
        ys_p = mesh_p.nodes[2]

        Ux_exact = [kovasznay_u(x, y) for x in xs_ux, y in ys_ux]
        Uy_exact = [kovasznay_v(x, y) for x in xs_uy, y in ys_uy]
        P_exact = [kovasznay_p(x, y) for x in xs_p, y in ys_p]

        x0_vec = zeros(total_dofs)
        x0_vec[1:nu_x] .= vec(Ux_exact)
        x0_vec[nu_x+1:2nu_x] .= vec(Ux_exact)
        x0_vec[2nu_x+1:2nu_x+nu_y] .= vec(Uy_exact)
        x0_vec[2nu_x+nu_y+1:2*(nu_x+nu_y)] .= vec(Uy_exact)
        x0_vec[2*(nu_x + nu_y) + 1:end] .= vec(P_exact)

        pressure_gauge = PinPressureGauge()
        interface_bc = Dirichlet(0.0)

        solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, interface_bc; x0=x0_vec)
        solve_NavierStokesMono_steady!(solver; tol=eps(), maxiter=8, relaxation=1.0, nlsolve_method=:picard)

        uωx = solver.x[1:nu_x]
        uωy = solver.x[2nu_x+1:2nu_x+nu_y]
        pω = solver.x[2*(nu_x + nu_y) + 1:end]

        Ux_num = reshape(uωx, Tuple(operator_ux.size))
        Uy_num = reshape(uωy, Tuple(operator_uy.size))
        P_num = reshape(pω, Tuple(operator_p.size))

        xs_ux_trim = trim_coords(xs_ux)
        ys_ux_trim = trim_coords(ys_ux)
        xs_uy_trim = trim_coords(xs_uy)
        ys_uy_trim = trim_coords(ys_uy)
        xs_p_trim = trim_coords(xs_p)
        ys_p_trim = trim_coords(ys_p)

        Ux_num_trim = trim_copy(Ux_num)
        Uy_num_trim = trim_copy(Uy_num)
        P_num_trim = trim_copy(P_num)

        Ux_exact_trim = trim_copy(Ux_exact)
        Uy_exact_trim = trim_copy(Uy_exact)
        P_exact_trim = trim_copy(P_exact)

        hx = Lx / nx
        hy = Ly / ny
        area_vel = hx * hy
        area_p = hx * hy

        ΔUx = Ux_num_trim .- Ux_exact_trim
        ΔUy = Uy_num_trim .- Uy_exact_trim
        ΔP = P_num_trim .- P_exact_trim
        ΔP .-= mean(ΔP)  # adjust pressure gauge

        err_u_L2 = sqrt(sum(abs2, ΔUx) * area_vel)
        err_v_L2 = sqrt(sum(abs2, ΔUy) * area_vel)
        err_p_L2 = sqrt(sum(abs2, ΔP) * area_p)

        err_u_Linf = maximum(abs, ΔUx)
        err_v_Linf = maximum(abs, ΔUy)
        err_p_Linf = maximum(abs, ΔP)

        h = max(hx, hy)

        push!(results, (; N, h, err_u_L2, err_u_Linf, err_v_L2, err_v_Linf, err_p_L2, err_p_Linf))

        println(@sprintf("%6d %7.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e",
                        N, h, err_u_L2, err_u_Linf, err_v_L2, err_v_Linf, err_p_L2, err_p_Linf))
    end

    return results
end

function convergence_rates(errors::Vector{Float64}, hs::Vector{Float64})
    rates = Float64[]
    for i in 1:length(errors)-1
        push!(rates, log(errors[i] / errors[i+1]) / log(hs[i] / hs[i+1]))
    end
    return rates
end

function write_results_csv(results; csv_path=nothing)
    hs = [r.h for r in results]
    u_L2_rates = convergence_rates([r.err_u_L2 for r in results], hs)
    v_L2_rates = convergence_rates([r.err_v_L2 for r in results], hs)
    p_L2_rates = convergence_rates([r.err_p_L2 for r in results], hs)
    u_Linf_rates = convergence_rates([r.err_u_Linf for r in results], hs)
    v_Linf_rates = convergence_rates([r.err_v_Linf for r in results], hs)
    p_Linf_rates = convergence_rates([r.err_p_Linf for r in results], hs)

    df = DataFrame(
        N = [r.N for r in results],
        h = hs,
        err_u_L2 = [r.err_u_L2 for r in results],
        err_u_Linf = [r.err_u_Linf for r in results],
        err_v_L2 = [r.err_v_L2 for r in results],
        err_v_Linf = [r.err_v_Linf for r in results],
        err_p_L2 = [r.err_p_L2 for r in results],
        err_p_Linf = [r.err_p_Linf for r in results],
        rate_u_L2 = vcat(NaN, u_L2_rates...),
        rate_u_Linf = vcat(NaN, u_Linf_rates...),
        rate_v_L2 = vcat(NaN, v_L2_rates...),
        rate_v_Linf = vcat(NaN, v_Linf_rates...),
        rate_p_L2 = vcat(NaN, p_L2_rates...),
        rate_p_Linf = vcat(NaN, p_Linf_rates...)
    )

    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "Kovasznay_2D_Mono.csv") : csv_path
    CSV.write(csv_file, df)
    return (csv_file=csv_file, table=df)
end

function main(; resolutions=[32, 64, 128, 256], csv_path=nothing)
    res = run_kovasznay(resolutions)
    csv_info = write_results_csv(res; csv_path=csv_path)
    return (results=res, csv_path=csv_info.csv_file, table=csv_info.table)
end

results = main()

@testset "Kovasznay 2D mono convergence" begin
    @test length(results.results) == 4
    @test results.results[1].h > results.results[end].h
    @test minimum([r.err_u_L2 for r in results.results]) < maximum([r.err_u_L2 for r in results.results])
    @test isfile(results.csv_path)
end

println("Kovasznay convergence complete. Results saved to ", results.csv_path)
