using Penguin
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using CSV
using DataFrames
using SpecialFunctions
using QuadGK
using Test

"""
Diphasic 2D heat diffusion benchmark (circular interface) that performs a
mesh-convergence study on mean interfacial flux and interfacial concentration.
Analytical reference values use the same integral expressions as the base
benchmark. Writes CSV data only (no plots).
"""
const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

const BASE_DL = 1.0
const DEFAULT_HE_DIFFUSIVITY_CASES = Dict(
    1e-3 => [1e0, 1e2, 1e4, 1e6],
    1e-2 => [1e0, 1e2, 1e4, 1e6],
    1e-1 => [1e0, 1e2, 1e4, 1e6],
    1e0  => [1e0, 1e2, 1e4, 1e6],
    1e1  => [1e0, 1e2, 1e4, 1e6],
    1e2  => [1e0, 1e2, 1e4, 1e6],
    1e3  => [1e0, 1e2, 1e4, 1e6]
)

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

Heat2Ph2DParams(; lx=10.0, ly=10.0, x0=0.0, y0=0.0, center=(5.0, 5.0),
                radius=1.0, Tend=0.1, Dg=1.0, Dl=1.0, He=1.0,
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

cg_prefactor(params::Heat2Ph2DParams) =
    (4 * params.cg0 * params.Dg * params.Dl^2 * params.He) / (π^2 * params.radius)

cl_prefactor(params::Heat2Ph2DParams) =
    (2 * params.cg0 * params.Dg * sqrt(params.Dl) * params.He) / π

function interfacial_concentrations(params::Heat2Ph2DParams; t=params.Tend, atol=1e-6, rtol=1e-6, Ufac=5.0)
    t > 0 || error("Need t > 0 for interfacial concentrations.")
    Umax = Ufac / sqrt(params.Dg * t)

    integrand_cg(u) = begin
        Φu = phi_val(u, params)
        Ψu = psi_val(u, params)
        denom = u^2 * (Φu^2 + Ψu^2)
        denom == 0.0 && return 0.0
        exp(-params.Dg * u^2 * t) * besselj0(u * params.radius) * besselj1(u * params.radius) / denom
    end

    integrand_cl(u) = begin
        Φu = phi_val(u, params)
        Ψu = psi_val(u, params)
        D = sqrt(params.Dg / params.Dl)
        denom = u * (Φu^2 + Ψu^2)
        denom == 0.0 && return 0.0
        contrib = besselj0(D * u * params.radius) * Φu - bessely0(D * u * params.radius) * Ψu
        exp(-params.Dg * u^2 * t) * besselj1(u * params.radius) * contrib / denom
    end

    Ig, _ = quadgk(integrand_cg, 0.0, Umax; atol=atol, rtol=rtol)
    Il, _ = quadgk(integrand_cl, 0.0, Umax; atol=atol, rtol=rtol)

    return (
        cgR = cg_prefactor(params) * Ig,
        clR = cl_prefactor(params) * Il,
        Umax = Umax
    )
end

function interfacial_flux_gas(params::Heat2Ph2DParams; t=params.Tend, atol=1e-6, rtol=1e-6)
    t > 0 || error("Need t > 0 for interfacial flux.")
    Ag = (4 * params.cg0 * params.Dg * params.Dl^2 * params.He) / (π^2 * params.radius)
    Umax = 5.0 / sqrt(params.Dg * t)

    integrand(u) = begin
        Φu = phi_val(u, params)
        Ψu = psi_val(u, params)
        denom = u * (Φu^2 + Ψu^2)
        denom == 0.0 && return 0.0
        exp(-params.Dg * u^2 * t) * besselj1(u * params.radius)^2 / denom
    end

    I, _ = quadgk(integrand, 0.0, Umax; atol=atol, rtol=rtol)
    return params.Dg * Ag * I
end

function interface_metrics(solver, operator1, operator2, capacity1, capacity2, params::Heat2Ph2DParams)
    n = prod(operator1.size)
    state = solver.states[end]
    Tω = @view state[1:n]
    Tγ = @view state[n+1:2n]
    Tω_c = @view state[2n+1:3n]
    Tγ_c = @view state[3n+1:4n]

    Q = operator1.H' * operator1.Wꜝ * (operator1.G * Tω + operator1.H * Tγ)
    Q_c = operator2.H' * operator2.Wꜝ * (operator2.G * Tω_c + operator2.H * Tγ_c)

    Γ = diag(capacity1.Γ)
    Γ_c = diag(capacity2.Γ)

    mask_Γ = .!isnan.(Γ)
    mask_Γ_c = .!isnan.(Γ_c)
    sum_Γ = sum(Γ[mask_Γ])
    sum_Γ_c = sum(Γ_c[mask_Γ_c])

    q_n_mean = sum(Q[mask_Γ]) / sum_Γ
    q_n_mean_c = sum(Q_c[mask_Γ_c]) / sum_Γ_c
    q_n_mean *= -params.Dg
    q_n_mean_c *= -params.Dl

    mask_iface = mask_Γ .& (Γ .> 0.0)
    mask_iface_c = mask_Γ_c .& (Γ_c .> 0.0)
    c_interface = sum(Tγ[mask_iface] .* Γ[mask_iface]) / sum(Γ[mask_iface])
    c_interface_c = sum(Tγ_c[mask_iface_c] .* Γ_c[mask_iface_c]) / sum(Γ_c[mask_iface_c])

    return (
        q_n_mean = q_n_mean,
        q_n_mean_c = q_n_mean_c,
        c_interface = c_interface,
        c_interface_c = c_interface_c
    )
end

error_metric(num, ana; relative=false) = begin
    rel = relative && ana != 0.0
    rel ? abs(num - ana) / abs(ana) : abs(num - ana)
end

function run_heat_2ph_2d_interface(
    nx_list::Vector{Int};
    params::Heat2Ph2DParams=Heat2Ph2DParams(),
    relative::Bool=false
)
    nx_vals = Int[]
    h_vals = Float64[]
    dt_vals = Float64[]

    qn_mean_gas = Float64[]
    qn_mean_liq = Float64[]
    qn_mean_ana_gas = Float64[]
    qn_mean_ana_liq = Float64[]
    qn_err_gas = Float64[]
    qn_err_liq = Float64[]
    qn_err_max = Float64[]

    c_interface_gas = Float64[]
    c_interface_liq = Float64[]
    c_interface_ana_gas = Float64[]
    c_interface_ana_liq = Float64[]
    c_interface_err_gas = Float64[]
    c_interface_err_liq = Float64[]
    c_interface_err_max = Float64[]

    for nx in nx_list
        ny = nx
        mesh = Penguin.Mesh((nx, ny), (params.lx, params.ly), (params.x0, params.y0))
        circle = (x, y, _=0) -> hypot(x - params.center[1], y - params.center[2]) - params.radius
        circle_c = (x, y, _=0) -> params.radius - hypot(x - params.center[1], y - params.center[2])

        capacity1 = Capacity(circle, mesh)
        capacity2 = Capacity(circle_c, mesh)
        operator1 = DiffusionOps(capacity1)
        operator2 = DiffusionOps(capacity2)

        bc_b = BorderConditions(Dict{Symbol,AbstractBoundary}())
        ic = InterfaceConditions(
            ScalarJump(params.He, 1.0, 0.0),
            FluxJump(params.Dg, params.Dl, 0.0)
        )

        f_zero = (x, y, z, t) -> 0.0
        D1_func = (x, y, z) -> params.Dg
        D2_func = (x, y, z) -> params.Dl

        phase1 = Phase(capacity1, operator1, f_zero, D1_func)
        phase2 = Phase(capacity2, operator2, f_zero, D2_func)

        ndofs = (nx + 1) * (ny + 1)
        u0 = vcat(ones(ndofs), ones(ndofs), zeros(ndofs), zeros(ndofs))
        Δt = 0.5 * (params.lx / nx)^2

        solver = DiffusionUnsteadyDiph(phase1, phase2, bc_b, ic, Δt, u0, "BE")
        solve_DiffusionUnsteadyDiph!(solver, phase1, phase2, Δt, params.Tend, bc_b, ic, "BE"; method=Base.:\)
        push!(solver.states, solver.x)

        metrics = interface_metrics(solver, operator1, operator2, capacity1, capacity2, params)
        ana_flux = interfacial_flux_gas(params)
        ana_flux_c = -ana_flux
        ana_conc = interfacial_concentrations(params)

        push!(nx_vals, nx)
        push!(h_vals, params.lx / nx)
        push!(dt_vals, Δt)

        push!(qn_mean_gas, metrics.q_n_mean)
        push!(qn_mean_liq, metrics.q_n_mean_c)
        push!(qn_mean_ana_gas, ana_flux)
        push!(qn_mean_ana_liq, ana_flux_c)

        err_gas = error_metric(metrics.q_n_mean, ana_flux; relative=relative)
        err_liq = error_metric(metrics.q_n_mean_c, ana_flux_c; relative=relative)
        push!(qn_err_gas, err_gas)
        push!(qn_err_liq, err_liq)
        push!(qn_err_max, max(err_gas, err_liq))

        push!(c_interface_gas, metrics.c_interface)
        push!(c_interface_liq, metrics.c_interface_c)
        push!(c_interface_ana_gas, ana_conc.cgR)
        push!(c_interface_ana_liq, ana_conc.clR)

        err_cg = error_metric(metrics.c_interface, ana_conc.cgR; relative=relative)
        err_cl = error_metric(metrics.c_interface_c, ana_conc.clR; relative=relative)
        push!(c_interface_err_gas, err_cg)
        push!(c_interface_err_liq, err_cl)
        push!(c_interface_err_max, max(err_cg, err_cl))
    end

    return (
        nx_vals = nx_vals,
        h_vals = h_vals,
        dt_vals = dt_vals,
        qn_mean_gas = qn_mean_gas,
        qn_mean_liq = qn_mean_liq,
        qn_mean_ana_gas = qn_mean_ana_gas,
        qn_mean_ana_liq = qn_mean_ana_liq,
        qn_err_gas = qn_err_gas,
        qn_err_liq = qn_err_liq,
        qn_err_max = qn_err_max,
        qn_orders = compute_orders(h_vals, qn_err_max, qn_err_max, qn_err_max),
        qn_pairwise_orders = compute_pairwise_orders(h_vals, qn_err_max),
        c_interface_gas = c_interface_gas,
        c_interface_liq = c_interface_liq,
        c_interface_ana_gas = c_interface_ana_gas,
        c_interface_ana_liq = c_interface_ana_liq,
        c_interface_err_gas = c_interface_err_gas,
        c_interface_err_liq = c_interface_err_liq,
        c_interface_err_max = c_interface_err_max,
        c_orders = compute_orders(h_vals, c_interface_err_max, c_interface_err_max, c_interface_err_max),
        c_pairwise_orders = compute_pairwise_orders(h_vals, c_interface_err_max),
        relative = relative
    )
end

function make_interface_convergence_dataframe(data)
    return DataFrame(
        nx = data.nx_vals,
        h = data.h_vals,
        dt = data.dt_vals,
        qn_mean_gas = data.qn_mean_gas,
        qn_mean_liq = data.qn_mean_liq,
        qn_mean_ana_gas = data.qn_mean_ana_gas,
        qn_mean_ana_liq = data.qn_mean_ana_liq,
        qn_err_gas = data.qn_err_gas,
        qn_err_liq = data.qn_err_liq,
        qn_err_max = data.qn_err_max,
        qn_pairwise_orders = data.qn_pairwise_orders,
        c_interface_gas = data.c_interface_gas,
        c_interface_liq = data.c_interface_liq,
        c_interface_ana_gas = data.c_interface_ana_gas,
        c_interface_ana_liq = data.c_interface_ana_liq,
        c_interface_err_gas = data.c_interface_err_gas,
        c_interface_err_liq = data.c_interface_err_liq,
        c_interface_err_max = data.c_interface_err_max,
        c_pairwise_orders = data.c_pairwise_orders
    )
end

function write_convergence_csv(method_name, data; csv_path=nothing)
    df = make_interface_convergence_dataframe(data)
    is_csv_file = !isnothing(csv_path) && endswith(lowercase(csv_path), ".csv")
    results_dir = isnothing(csv_path) ?
        joinpath(BENCH_ROOT, "results", "scalar", "diphasic", "extreme_regimes", "interface") :
        (is_csv_file ? dirname(csv_path) : csv_path)
    mkpath(results_dir)
    csv_out = (isnothing(csv_path) || !is_csv_file) ?
        joinpath(results_dir, "$(method_name)_Convergence.csv") :
        csv_path
    CSV.write(csv_out, df)
    return (csv_path = csv_out, table = df)
end

function _rebuild_params(params::Heat2Ph2DParams; ratio=1.0, Dl=BASE_DL, He=params.He)
    return Heat2Ph2DParams(
        params.lx, params.ly, params.x0, params.y0, params.center, params.radius,
        params.Tend, ratio * Dl, Dl, He, params.cg0, params.cl0
    )
end

function main(;
    csv_path=nothing,
    csv_dir=csv_path,
    nx_list=nothing,
    case_dict=nothing,
    Dl=BASE_DL,
    params::Heat2Ph2DParams=Heat2Ph2DParams(),
    relative=false
)
    nx_vals = isnothing(nx_list) ? [8, 16, 32, 64, 128, 256] : nx_list
    cases = isnothing(case_dict) ? DEFAULT_HE_DIFFUSIVITY_CASES : case_dict

    if length(cases) > 1 && !isnothing(csv_dir) && endswith(lowercase(csv_dir), ".csv")
        error("Provide a directory path (or nothing) for csv_dir when running multiple cases.")
    end

    results = NamedTuple[]

    for (He, ratios) in cases
        for ratio in ratios
            params_case = _rebuild_params(params; ratio=ratio, Dl=Dl, He=He)
            data = run_heat_2ph_2d_interface(nx_vals; params=params_case, relative=relative)
            method_name = "Heat_2ph_2D_interface_ratio$(ratio)_He$(params_case.He)_Dl$(params_case.Dl)"
            csv_info = write_convergence_csv(method_name, data; csv_path=csv_dir)
            push!(results, (
                ratio = ratio,
                He = params_case.He,
                params = params_case,
                data = data,
                csv_path = csv_info.csv_path,
                table = csv_info.table
            ))
        end
    end

    return results
end

results = main()

@testset "Diphasic Heat 2D interface convergence" begin
    expected = sum(length(v) for v in values(DEFAULT_HE_DIFFUSIVITY_CASES))
    @test length(results) == expected
    for res in results
        @test length(res.data.h_vals) == length(res.data.qn_err_max)
        @test length(res.data.h_vals) == length(res.data.c_interface_err_max)
        @test res.data.h_vals[1] > res.data.h_vals[end]
        @test minimum(res.data.qn_err_max) <= maximum(res.data.qn_err_max)
        @test minimum(res.data.c_interface_err_max) <= maximum(res.data.c_interface_err_max)
        @test isfile(res.csv_path)
    end
end
