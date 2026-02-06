using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using SpecialFunctions
using QuadGK
using Test
using Printf

"""
Diphasic 2D heat diffusion benchmark reproduced from `benchmark/Heat_2ph_2D.jl`.
Performs a mesh-convergence study with a circular interface and writes CSV data
only (no plots). Custom parameter set with Δt = Tend / N.
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

struct Heat2Ph2DParams
    lx::Float64
    ly::Float64
    x0::Float64
    y0::Float64
    center::Tuple{Float64,Float64}
    radius::Float64
    Tend::Float64
    Dg::Float64
    Dl::Float64
    He::Float64
    cg0::Float64
    cl0::Float64
end

Heat2Ph2DParams(; lx=16.0, ly=16.0, x0=0.0, y0=0.0, center=(8.0+0.125/pi, 8.0+0.125*sqrt(3/2)),
                radius=1.0, Tend=3.0, Dg=17.416, Dl=0.17416, He=1/30.0,
                cg0=1.0, cl0=0.0) =
    Heat2Ph2DParams(lx, ly, x0, y0, center, radius, Tend, Dg, Dl, He, cg0, cl0)

phi_val(u, params::Heat2Ph2DParams) = begin
    D = sqrt(params.Dg / params.Dl)
    term1 = params.Dg * sqrt(params.Dl) * besselj1(u * params.radius) * bessely0(D * u * params.radius)
    term2 = params.He * params.Dl * sqrt(params.Dg) * besselj0(u * params.radius) * bessely1(D * u * params.radius)
    term1 - term2
end

psi_val(u, params::Heat2Ph2DParams) = begin
    D = sqrt(params.Dg / params.Dl)
    term1 = params.Dg * sqrt(params.Dl) * besselj1(u * params.radius) * besselj0(D * u * params.radius)
    term2 = params.He * params.Dl * sqrt(params.Dg) * besselj0(u * params.radius) * besselj1(D * u * params.radius)
    term1 - term2
end

function cg_integrand(u, x, y, params::Heat2Ph2DParams)
    r = hypot(x - params.center[1], y - params.center[2])
    Φu = phi_val(u, params)
    Ψu = psi_val(u, params)
    denom = u^2 * (Φu^2 + Ψu^2)
    numer = exp(-params.Dg * u^2 * params.Tend) * besselj0(u * r) * besselj1(u * params.radius)
    iszero(denom) ? 0.0 : numer / denom
end

function cl_integrand(u, x, y, params::Heat2Ph2DParams)
    r = hypot(x - params.center[1], y - params.center[2])
    Φu = phi_val(u, params)
    Ψu = psi_val(u, params)
    D = sqrt(params.Dg / params.Dl)
    denom = u * (Φu^2 + Ψu^2)
    contrib = besselj0(D * u * r) * Φu - bessely0(D * u * r) * Ψu
    numer = exp(-params.Dg * u^2 * params.Tend) * besselj1(u * params.radius) * contrib
    iszero(denom) ? 0.0 : numer / denom
end

function heat2d_phase1_solution(params::Heat2Ph2DParams)
    prefactor = (4 * params.cg0 * params.Dg * params.Dl^2 * params.He) / (π^2 * params.radius)
    Umax = 5.0 / sqrt(params.Dg * params.Tend)
    return (x, y) -> begin
        r = hypot(x - params.center[1], y - params.center[2])
        r >= params.radius && return 0.0
        val, _ = quadgk(u -> cg_integrand(u, x, y, params), 0, Umax; atol=1e-6, rtol=1e-6)
        prefactor * val
    end
end

function heat2d_phase2_solution(params::Heat2Ph2DParams)
    prefactor = (2 * params.cg0 * params.Dg * sqrt(params.Dl) * params.He) / π
    Umax = 5.0 / sqrt(params.Dg * params.Tend)
    return (x, y) -> begin
        r = hypot(x - params.center[1], y - params.center[2])
        r < params.radius && return 0.0
        val, _ = quadgk(u -> cl_integrand(u, x, y, params), 0, Umax; atol=1e-6, rtol=1e-6)
        prefactor * val
    end
end

function save_cell_data(nx, solver, capacity1, capacity2, u1_exact, u2_exact, params::Heat2Ph2DParams, output_dir)
    """Save cell data for a given mesh resolution"""
    # Extract solution vectors - inside (phase1) and outside (phase2) solutions
    ndofs = length(capacity1.C_ω)
    state = solver.states[end]
    cg_inside = state[1:ndofs]                 # phase 1 bulk values
    cl_outside = state[2*ndofs+1:3*ndofs]       # phase 2 bulk values

    # Get cell centroids from both capacities
    centroids1 = capacity1.C_ω
    centroids2 = capacity2.C_ω

    # Accumulate weighted L2 (mean-square) errors using uniform ΔxΔy
    dx = capacity1.mesh.nodes[1][2] - capacity1.mesh.nodes[1][1]
    dy = capacity1.mesh.nodes[2][2] - capacity1.mesh.nodes[2][1]
    dA = dx * dy
    sum_w_inside = 0.0
    sum_werr2_inside = 0.0
    sum_w_outside = 0.0
    sum_werr2_outside = 0.0

    # Save inside cylinder data (drop only the last row)
    open(joinpath(output_dir, "mesh_nx$(nx)_inside_cell_data.txt"), "w") do io
        write(io, "x_pos\ty_pos\tradius\tcg_num\tcg_exa\n")

        for cell_idx in 1:(length(centroids1) - 1)
            centroid = centroids1[cell_idx]
            x_center = centroid[1]
            y_center = centroid[2]

            r = hypot(x_center - params.center[1], y_center - params.center[2])
            cg_val = cg_inside[cell_idx]
            cg_exact_val = u1_exact(x_center, y_center)

            @printf(io, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                    x_center, y_center, r, cg_val, cg_exact_val)

            err = cg_exact_val - cg_val
            sum_w_inside += dA
            sum_werr2_inside += dA * err * err
        end
    end

    # Save outside cylinder data (drop only the last row)
    open(joinpath(output_dir, "mesh_nx$(nx)_outside_cell_data.txt"), "w") do io
        write(io, "x_pos\ty_pos\tradius\tcl_num\tcl_exa\n")

        for cell_idx in 1:(length(centroids2) - 1)
            centroid = centroids2[cell_idx]
            x_center = centroid[1]
            y_center = centroid[2]

            r = hypot(x_center - params.center[1], y_center - params.center[2])
            cl_val = cl_outside[cell_idx]
            cl_exact_val = u2_exact(x_center, y_center)

            @printf(io, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                    x_center, y_center, r, cl_val, cl_exact_val)

            err = cl_exact_val - cl_val
            sum_w_outside += dA
            sum_werr2_outside += dA * err * err
        end
    end

    l2_inside = sum_w_inside > 0 ? sqrt(sum_werr2_inside / sum_w_inside) : NaN
    l2_outside = sum_w_outside > 0 ? sqrt(sum_werr2_outside / sum_w_outside) : NaN
    open(joinpath(output_dir, "mesh_nx$(nx)_l2_norms.txt"), "w") do io
        write(io, "phase\tl2_rms\n")
        @printf(io, "inside\t%.6e\n", l2_inside)
        @printf(io, "outside\t%.6e\n", l2_outside)
    end
end

function save_interface_cell_data(nx, solver, capacity1, capacity2, operator1, operator2, params::Heat2Ph2DParams, output_dir)
    """Save interface flux (Q, Q_c) and interfacial concentrations (Tgamma, Tgamma_c)"""
    n = prod(operator1.size)
    state = solver.states[end]
    Tω = @view state[1:n]
    Tgamma = @view state[n+1:2n]
    Tω_c = @view state[2n+1:3n]
    Tgamma_c = @view state[3n+1:4n]

    Q = params.Dg * operator1.H' * operator1.Wꜝ * (operator1.G * Tω + operator1.H * Tgamma)
    Q_c = params.Dl * operator2.H' * operator2.Wꜝ * (operator2.G * Tω_c + operator2.H * Tgamma_c)

    Γ = diag(capacity1.Γ)
    Γ_c = diag(capacity2.Γ)
    mask_iface = .!isnan.(Γ) .& (Γ .> 0.0)
    mask_iface_c = .!isnan.(Γ_c) .& (Γ_c .> 0.0)

    centroids1 = capacity1.C_γ
    centroids2 = capacity2.C_γ

    # use the word inside
    open(joinpath(output_dir, "mesh_nx$(nx)_interface_inside_data.txt"), "w") do io
        write(io, "x_pos\ty_pos\tradius\tGamma\tTgamma\tQ\n")
        for idx in eachindex(Γ)
            mask_iface[idx] || continue
            centroid = centroids1[idx]
            x_center = centroid[1]
            y_center = centroid[2]
            r = hypot(x_center - params.center[1], y_center - params.center[2])
            @printf(io, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                    x_center, y_center, r, Γ[idx], Tgamma[idx], Q[idx])
        end
    end

    open(joinpath(output_dir, "mesh_nx$(nx)_interface_outside_data.txt"), "w") do io
        write(io, "x_pos\ty_pos\tradius\tGamma\tTgamma_c\tQ_c\n")
        for idx in eachindex(Γ_c)
            mask_iface_c[idx] || continue
            centroid = centroids2[idx]
            x_center = centroid[1]
            y_center = centroid[2]
            r = hypot(x_center - params.center[1], y_center - params.center[2])
            @printf(io, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                    x_center, y_center, r, Γ_c[idx], Tgamma_c[idx], Q_c[idx])
        end
    end
end

function run_heat_2ph_2d(
    nx_list::Vector{Int};
    params::Heat2Ph2DParams=Heat2Ph2DParams(),
    norm::Real=2,
    relative::Bool=false,
    save_cell_data_flag::Bool=true
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

    u1_exact = heat2d_phase1_solution(params)
    u2_exact = heat2d_phase2_solution(params)
    
    # Create output directory for cell data
    cell_data_dir = joinpath(BENCH_ROOT, "results", "scalar", "diphasic", "extreme_regimes", "cell_data_shifted_large_He$(params.He)_Ratio$(params.Dg/params.Dl)_int")
    if save_cell_data_flag
        mkpath(cell_data_dir)
    end

    for nx in nx_list
        ny = nx
        mesh = Penguin.Mesh((nx, ny), (params.lx, params.ly), (params.x0, params.y0))
        circle = (x, y, _=0) -> sqrt((x - params.center[1])^2 + (y - params.center[2])^2) - params.radius
        circle_c = (x, y, _=0) -> params.radius - sqrt((x - params.center[1])^2 + (y - params.center[2])^2)
        capacity1 = Capacity(circle, mesh; compute_centroids=true)
        capacity2 = Capacity(circle_c, mesh; compute_centroids=true)
        operator1 = DiffusionOps(capacity1)
        operator2 = DiffusionOps(capacity2)

        bc_b = BorderConditions(Dict{Symbol,AbstractBoundary}())
        ic = InterfaceConditions(
            ScalarJump(1.0, 1/params.He, 0.0),
            FluxJump(params.Dg, params.Dl, 0.0)
        )

        f_zero = (x, y, z, t) -> 0.0
        D1_func = (x, y, z) -> params.Dg
        D2_func = (x, y, z) -> params.Dl

        phase1 = Phase(capacity1, operator1, f_zero, D1_func)
        phase2 = Phase(capacity2, operator2, f_zero, D2_func)

        ndofs = (nx + 1) * (ny + 1)
        u0 = vcat(ones(ndofs), ones(ndofs), zeros(ndofs), zeros(ndofs))
        Δt = params.Tend / (2nx)
        push!(dt_vals, Δt)

        solver = DiffusionUnsteadyDiph(phase1, phase2, bc_b, ic, Δt, u0, "CN")
        solve_DiffusionUnsteadyDiph!(solver, phase1, phase2, Δt, params.Tend, bc_b, ic, "CN"; method=Base.:\)
        push!(solver.states, solver.x)

        _, _, global_errs, full_errs, cut_errs, empty_errs =
            check_convergence_diph(u1_exact, u2_exact, solver, capacity1, capacity2, norm, relative)

        push!(h_vals, params.lx / nx)
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
        
        # Save cell data for this mesh
        if save_cell_data_flag
            save_cell_data(nx, solver, capacity1, capacity2, u1_exact, u2_exact, params, cell_data_dir)
            save_interface_cell_data(nx, solver, capacity1, capacity2, operator1, operator2, params, cell_data_dir)
        end
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
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic", "extreme_regimes") :
        dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ?
        joinpath(results_dir, "$(method_name)_Convergence.csv") :
        csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function main(; csv_path=nothing, nx_list=nothing, params::Heat2Ph2DParams=Heat2Ph2DParams())
    nx_vals = isnothing(nx_list) ? [16, 32, 64, 128, 256, 512] : nx_list
    data = run_heat_2ph_2d(nx_vals; params=params)
    csv_info = write_convergence_csv("Heat_2ph_2D_phase_custom_He$(params.He)_Ratio$(params.Dg/params.Dl)", data; csv_path=csv_path)
    return (data = data, csv_path = csv_info.csv_path, table = csv_info.table)
end

results = main()

@testset "Diphasic Heat 2D convergence (phase custom regime)" begin
    orders = results.data.orders
    #@test !isnan(orders.all)
    @test length(results.data.h_vals) == length(results.data.err_vals)
    @test results.data.h_vals[1] > results.data.h_vals[end]
    @test minimum(results.data.err_vals) < maximum(results.data.err_vals)
    @test isfile(results.csv_path)
end
