using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using QuadGK
using CSV
using Test

"""
Diphasic 3D diffusion benchmark with an analytical radial solution.

Domain: [0, 4]^3. Interface: sphere centered at (2, 2, 2) with radius R0 = 1.
Diffusivities: Dg = 1.0 (interior), Dl = 4.0 (exterior). Henry coefficient
alpha = 1, so value and flux are continuous across the interface.

Exact transient solution (bubble model):
  u_g(t, r) = -2 c0 / (pi * r) * int_0^{Xmax} imag((sinh(lambda_g r) / zeta(x)) * exp(-x * t)) / x dx
  u_l(t, r) =  c0 / (pi * r) * int_0^{Xmax} imag((xi(x) * exp(-lambda_l r) / zeta(x)) * exp(-x * t)) / x dx
where lambda_g = im * sqrt(x / Dg), lambda_l = im * sqrt(x / Dl),
      xi(x)    = 2 * Dg * (lambda_g * R0 * cosh(lambda_g * R0) - sinh(lambda_g * R0)) /
                 (Dl * exp(-lambda_l * R0) * (1 + lambda_l * R0)),
      zeta(x)  = xi(x) * exp(-lambda_l * R0) / (alpha * R0) + 2 * sinh(lambda_g * R0) / R0.

Boundary conditions: Dirichlet set to the exact solution on the box.
Source term: zero in both phases.
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

const Dg = 1.0
const Dl = 1.0
const R0 = 1.0
const alpha = 1.0        # Henry coefficient α
const c0 = 1.0           # Reference temperature T0
const X_MAX = 1000.0
const CENTER = (2.0, 2.0, 2.0)

# Parameters derived from the analytical solution
σ = sqrt(Dl / Dg)
Q = (Dl / Dg) * σ
L = (Dl - Dg) / Dg

const t1_cache = Dict{Tuple{Float64,Float64},Float64}()
const t2_cache = Dict{Tuple{Float64,Float64},Float64}()

denom(u) = (u * cos(u) + L * sin(u))^2 + (Q * u * sin(u))^2

function integrand_core(u, r, t)
    # Eq. (7) inner phase
    num = (sin(u) - u * cos(u)) * sin(u * r / R0)
    expo = exp(-Dg * u^2 * t / R0^2)
    return num * expo / denom(u)
end

function F_shell(u, r)
    arg = u * (r - R0) / (σ * R0)
    return (u * cos(u) + L * sin(u)) * sin(arg) + Q * u * sin(u) * cos(arg)
end

function integrand_shell(u, r, t)
    # Eq. (8) outer phase; include σ/u factor from du/u′ with u′=u/σ
    num = (sin(u) - u * cos(u)) * F_shell(u, r)
    expo = exp(-Dg * u^2 * t / R0^2)
    u_safe = max(u, 1e-12)
    return num * expo * σ / (u_safe * denom(u))
end

function T_core_exact(t, r)
    r_eff = max(r, 1e-8)
    key = (round(t, sigdigits=6), round(r_eff, sigdigits=6))
    return get!(t1_cache, key) do
        val, _ = quadgk(u -> integrand_core(u, r_eff, t), 0.0, X_MAX; atol=1e-6, rtol=1e-6, maxevals=10_000)
        (2 * Q * alpha / (pi * r_eff)) * val * c0
    end
end

function T_shell_exact(t, r)
    r_eff = max(r, 1e-8)
    key = (round(t, sigdigits=6), round(r_eff, sigdigits=6))
    return get!(t2_cache, key) do
        val, _ = quadgk(u -> integrand_shell(u, r_eff, t), 0.0, X_MAX; atol=1e-6, rtol=1e-6, maxevals=10_000)
        (2 * alpha / (pi * r_eff)) * val * c0
    end
end

exact_field(t, x, y, z; center=CENTER) = begin
    r = sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2)
    r < R0 ? T_core_exact(t, r) : T_shell_exact(t, r)
end

f_zero(x, y, z, t) = 0.0
D_interior(x, y, z) = Dg
D_exterior(x, y, z) = Dl

function build_bubble_components(mesh; center=CENTER, radius=R0)
    body_out = (x, y, z, _=0.0) -> sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2) - radius
    body_in  = (x, y, z, _=0.0) -> radius - sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2)

    capacity1 = Capacity(body_in, mesh; method="VOFI")
    capacity2 = Capacity(body_out, mesh; method="VOFI")
    operator1 = DiffusionOps(capacity1)
    operator2 = DiffusionOps(capacity2)

    bc_func = function (x, y, z, t=0.0)
        tval = isnothing(t) ? 0.0 : t
        return exact_field(tval, x, y, z; center=center)
    end

    bc_b = BorderConditions(Dict(
    ))

    ic = InterfaceConditions(
        ScalarJump(1.0, 1.0/alpha, 0.0),
        FluxJump(Dg, Dl, 0.0)
    )

    phase1 = Phase(capacity1, operator1, f_zero, D_interior)
    phase2 = Phase(capacity2, operator2, f_zero, D_exterior)

    return capacity1, capacity2, operator1, operator2, phase1, phase2, bc_b, ic
end

function run_bubble_3d(
    nx_list::Vector{Int};
    lx::Float64 = 8.0,
    ly::Float64 = 8.0,
    lz::Float64 = 8.0,
    center::Tuple{Float64,Float64,Float64} = CENTER,
    radius::Float64 = R0,
    Tend::Float64 = 0.1,
    norm::Real = 2,
    relative::Bool = false
)
    h_vals = Float64[]
    dt_vals = Float64[]
    err_vals = Float64[]
    err_full_vals = Float64[]
    err_cut_vals = Float64[]
    err_empty_vals = Float64[]
    phase1_all_errs = Float64[]
    phase1_full_errs = Float64[]
    phase1_cut_errs = Float64[]
    phase1_empty_errs = Float64[]
    phase2_all_errs = Float64[]
    phase2_full_errs = Float64[]
    phase2_cut_errs = Float64[]
    phase2_empty_errs = Float64[]
    inside_cells = Int[]
    inside_cells_by_dim = Vector{Vector{Int}}()
    inside_cells_phase1 = Int[]
    inside_cells_phase2 = Int[]

    for nx in nx_list
        ny = nx; nz = nx
        mesh = Penguin.Mesh((nx, ny, nz), (lx, ly, lz), (0.0, 0.0, 0.0))

        capacity1, capacity2, operator1, operator2, phase1, phase2, bc_b, ic =
            build_bubble_components(mesh; center=center, radius=radius)
        ndofs1 = length(capacity1.C_ω)
        ndofs2 = length(capacity2.C_ω)
        @assert ndofs1 == ndofs2 "Expected matching DOF counts for both phases"

        # Initial condition: exact field at t=0
        u0_phase1 = [exact_field(0.0, c[1], c[2], c[3]; center=center) for c in capacity1.C_ω]
        u0_phase2 = [exact_field(0.0, c[1], c[2], c[3]; center=center) for c in capacity2.C_ω]
        u0 = vcat(u0_phase1, u0_phase1, u0_phase2, u0_phase2)

        Δt = 0.5 * (minimum((lx / nx, ly / ny, lz / nz))^2) / max(Dg, Dl)
        push!(dt_vals, Δt)

        solver = DiffusionUnsteadyDiph(phase1, phase2, bc_b, ic, Δt, u0, "BE")
        solve_DiffusionUnsteadyDiph!(solver, phase1, phase2, Δt, Tend, bc_b, ic, "BE"; method=Base.:\)
        push!(solver.states, solver.x)

        u1_exact = (x, y, z) -> exact_field(Tend, x, y, z; center=center)
        u2_exact = u1_exact

        _, _, global_errs, full_errs, cut_errs, empty_errs =
            check_convergence_diphh(u1_exact, u2_exact, solver, capacity1, capacity2, norm, relative)

        push!(h_vals, minimum((lx / nx, ly / ny, lz / nz)))
        push!(err_vals, global_errs[3])
        push!(err_full_vals, full_errs[3])
        push!(err_cut_vals, cut_errs[3])
        push!(err_empty_vals, empty_errs[3])
        push!(phase1_all_errs, global_errs[1])
        push!(phase2_all_errs, global_errs[2])
        push!(phase1_full_errs, full_errs[1])
        push!(phase2_full_errs, full_errs[2])
        push!(phase1_cut_errs, cut_errs[1])
        push!(phase2_cut_errs, cut_errs[2])
        push!(phase1_empty_errs, empty_errs[1])
        push!(phase2_empty_errs, empty_errs[2])

        inside1 = count_inside_cells(capacity1)
        inside2 = count_inside_cells(capacity2)
        push!(inside_cells, inside1 + inside2)
        push!(inside_cells_by_dim, [inside1, inside2])
        push!(inside_cells_phase1, inside1)
        push!(inside_cells_phase2, inside2)
    end

    return (
        h_vals = h_vals,
        dt_vals = dt_vals,
        err_vals = err_vals,
        err_full_vals = err_full_vals,
        err_cut_vals = err_cut_vals,
        err_empty_vals = err_empty_vals,
        phase1_all_errs = phase1_all_errs,
        phase1_full_errs = phase1_full_errs,
        phase1_cut_errs = phase1_cut_errs,
        phase1_empty_errs = phase1_empty_errs,
        phase2_all_errs = phase2_all_errs,
        phase2_full_errs = phase2_full_errs,
        phase2_cut_errs = phase2_cut_errs,
        phase2_empty_errs = phase2_empty_errs,
        inside_cells = inside_cells,
        inside_cells_by_dim = inside_cells_by_dim,
        inside_cells_phase1 = inside_cells_phase1,
        inside_cells_phase2 = inside_cells_phase2,
        orders = compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals),
        norm = norm
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_diphasic_convergence_dataframe(method_name, data)
    results_dir = isnothing(csv_path) ?
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic") :
        dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ?
        joinpath(results_dir, "$(method_name)_Convergence.csv") :
        csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing)
    nx_vals = isnothing(nx_list) ? [4, 8, 16, 31, 42, 53] : nx_list
    data = run_bubble_3d(nx_vals)
    csv_info = write_convergence_csv("Diph_3D_Bubble", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Diphasic 3D bubble convergence" begin
    orders = results.data.orders
    @test !isnan(orders.all)
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end
