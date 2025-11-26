using Penguin
using Statistics
using LinearAlgebra
using CSV
using DataFrames
using Printf
using Test

"""
Oscillatory channel flow (Womersley-type) in 2D using the unsteady
Navier–Stokes mono solver. Produces CSV outputs only (no plots).
"""

const h = 1.0
const Lx = 2 * h
const nx = 64
const ny = 64
const ρ = 1.0
const ν = 0.05
const μ = ρ * ν

const y_wall_bot = -h
const y_wall_top = h

mesh_p  = Penguin.Mesh((nx, ny), (Lx, 2 * h), (0.0, -h))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (Lx, 2 * h), (-0.5 * dx, -h))
mesh_uy = Penguin.Mesh((nx, ny), (Lx, 2 * h), (0.0, -h - 0.5 * dy))

body = (x, y, _=0.0) -> begin
    if y < y_wall_bot
        return y_wall_bot - y
    elseif y > y_wall_top
        return y - y_wall_top
    else
        return -min(y - y_wall_bot, y_wall_top - y)
    end
end

capacity_ux = Capacity(body, mesh_ux; compute_centroids=false)
capacity_uy = Capacity(body, mesh_uy; compute_centroids=false)
capacity_p  = Capacity(body, mesh_p;  compute_centroids=false)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

nu_x = prod(operator_ux.size)
nu_y = prod(operator_uy.size)
np = prod(operator_p.size)
total_dofs = 2 * (nu_x + nu_y) + np

xs_ux = mesh_ux.nodes[1]
ys_ux = mesh_ux.nodes[2]
xs_uy = mesh_uy.nodes[1]
ys_uy = mesh_uy.nodes[2]
xs_p  = mesh_p.nodes[1]
ys_p  = mesh_p.nodes[2]

trim_coords(v) = length(v) > 1 ? v[1:end-1] : copy(v)
xs_ux_trim = trim_coords(xs_ux)
ys_ux_trim = trim_coords(ys_ux)
xs_uy_trim = trim_coords(xs_uy)
ys_uy_trim = trim_coords(ys_uy)

###########
# Utilities
###########
function complex_profile(F0, ω, ν_local, h_local, y)
    λ = sqrt(im * ω / ν_local)
    F0 / (im * ω) * (1 - cosh(λ * y) / cosh(λ * h_local))
end

function complex_bulk(F0, ω, ν_local, h_local)
    λ = sqrt(im * ω / ν_local)
    F0 / (im * ω) * (1 - tanh(λ * h_local) / (λ * h_local))
end

function complex_wall_shear(F0, ω, ν_local, h_local, μ_local)
    λ = sqrt(im * ω / ν_local)
    Uprime = -(F0 / (im * ω)) * λ * tanh(λ * h_local)
    μ_local * Uprime
end

function forcing_from_bulk(target_bulk, ω, ν_local, h_local)
    λ = sqrt(im * ω / ν_local)
    im * ω * target_bulk / (1 - tanh(λ * h_local) / (λ * h_local))
end

function harmonic_fit(t, s, ω)
    A = hcat(cos.(ω .* t), sin.(ω .* t))
    coeffs = A \ s
    coeffs[1] - im * coeffs[2]
end

wrap_phase_deg(θ) = mod(θ + 180, 360) - 180

function trapz(x, y)
    n = length(x)
    n == length(y) || error("trapz dimension mismatch")
    n == 1 && return 0.0
    acc = 0.0
    for i in 2:n
        acc += 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    end
    acc
end

analytic_velocity(y, t, F0, ω, ν_local, h_local) =
    real(complex_profile(F0, ω, ν_local, h_local, y) * exp(im * ω * t))

function analyze_case(alpha; target_bulk=1.0, periods=1, steps_per_period=240)
    ω = (alpha^2 * ν) / h^2
    T = 2π / ω
    Δt = T / steps_per_period
    total_steps = periods * steps_per_period
    T_end = periods * T

    F0 = forcing_from_bulk(target_bulk, ω, ν, h)
    U_bulk_exact = complex_bulk(F0, ω, ν, h)
    U_center_exact = complex_profile(F0, ω, ν, h, 0.0)
    τ_exact = complex_wall_shear(F0, ω, ν, h, μ)

    inlet_velocity = (x, y, t=0.0) -> analytic_velocity(y, t, F0, ω, ν, h)

    bc_ux = BorderConditions(Dict(
        :left   => Dirichlet(inlet_velocity),
        :right  => Outflow(),
        :bottom => Dirichlet((x, y, t=0.0) -> 0.0),
        :top    => Dirichlet((x, y, t=0.0) -> 0.0)
    ))

    bc_uy = BorderConditions(Dict(
        :left   => Dirichlet((x, y, t=0.0) -> 0.0),
        :right  => Outflow(),
        :bottom => Dirichlet((x, y, t=0.0) -> 0.0),
        :top    => Dirichlet((x, y, t=0.0) -> 0.0)
    ))

    f_force = (x, y, z=0.0, t=0.0) -> real(F0 * exp(im * ω * t))
    fluid = Fluid((mesh_ux, mesh_uy),
                  (capacity_ux, capacity_uy),
                  (operator_ux, operator_uy),
                  mesh_p,
                  capacity_p,
                  operator_p,
                  μ, ρ, f_force, (x, y, z=0.0, t=0.0)->0.0)

    pressure_gauge = PinPressureGauge()
    interface_bc = Dirichlet(0.0)

    solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, interface_bc; x0=zeros(total_dofs))
    times, histories = solve_NavierStokesMono_unsteady!(solver; Δt=Δt, T_end=T_end, scheme=:CN)

    xs_len = length(xs_ux_trim)
    ys_len = length(ys_ux_trim)
    center_idx = clamp(argmin(abs.(ys_ux_trim)), 1, ys_len)

    total_len = length(times)
    center_series = zeros(Float64, total_len)
    bulk_series = zeros(Float64, total_len)
    shear_series = zeros(Float64, total_len)
    force_series = zeros(Float64, total_len)
    power_series = zeros(Float64, total_len)
    diss_series = zeros(Float64, total_len)

    profile_series = zeros(Float64, ys_len, steps_per_period)
    profile_counter = 0

    for (step, hist) in enumerate(histories)
        ux_hist = hist[1:nu_x]
        Ux = reshape(ux_hist, (length(xs_ux), length(ys_ux)))
        Ux_trim = Ux[1:length(xs_ux_trim), 1:length(ys_ux_trim)]
        avg_profile = vec(mean(Ux_trim, dims=1))

        center_series[step] = avg_profile[center_idx]
        bulk_series[step] = trapz(ys_ux_trim, avg_profile) / (2 * h)

        du_dy = similar(avg_profile)
        du_dy[1] = (avg_profile[2] - avg_profile[1]) / (ys_ux_trim[2] - ys_ux_trim[1])
        du_dy[end] = (avg_profile[end] - avg_profile[end-1]) / (ys_ux_trim[end] - ys_ux_trim[end-1])
        for j in 2:length(avg_profile)-1
            du_dy[j] = (avg_profile[j+1] - avg_profile[j-1]) / (ys_ux_trim[j+1] - ys_ux_trim[j-1])
        end
        shear_series[step] = μ * du_dy[end]

        f_val = real(F0 * exp(im * ω * times[step]))
        force_series[step] = f_val
        power_series[step] = ρ * f_val * trapz(ys_ux_trim, avg_profile) * Lx
        diss_series[step] = μ * trapz(ys_ux_trim, du_dy .^ 2) * Lx

        if step > total_len - steps_per_period
            profile_counter += 1
            profile_series[:, profile_counter] = avg_profile
        end
    end

    idx_start = total_len - steps_per_period + 1
    t_window = times[idx_start:end]
    center_window = center_series[idx_start:end]
    bulk_window = bulk_series[idx_start:end]
    shear_window = shear_series[idx_start:end]
    force_window = force_series[idx_start:end]
    power_window = power_series[idx_start:end]
    diss_window = diss_series[idx_start:end]

    center_amp_num = harmonic_fit(t_window, center_window, ω)
    bulk_amp_num = harmonic_fit(t_window, bulk_window, ω)
    shear_amp_num = harmonic_fit(t_window, shear_window, ω)

    amp_err_center = abs(abs(center_amp_num) - abs(U_center_exact)) / max(abs(U_center_exact), eps())
    amp_err_bulk = abs(abs(bulk_amp_num) - abs(U_bulk_exact)) / max(abs(U_bulk_exact), eps())
    amp_err_shear = abs(abs(shear_amp_num) - abs(τ_exact)) / max(abs(τ_exact), eps())

    phase_err_center = wrap_phase_deg(rad2deg(angle(center_amp_num / U_center_exact)))
    phase_err_bulk = wrap_phase_deg(rad2deg(angle(bulk_amp_num / U_bulk_exact)))
    phase_err_shear = wrap_phase_deg(rad2deg(angle(shear_amp_num / τ_exact)))

    profile_amplitudes = Vector{ComplexF64}(undef, ys_len)
    for j in 1:ys_len
        profile_amplitudes[j] = harmonic_fit(t_window, view(profile_series, j, :), ω)
    end
    profile_exact = [complex_profile(F0, ω, ν, h, y) for y in ys_ux_trim]

    profile_tag = replace(@sprintf("alpha_%g", alpha), "." => "p")
    profile_df = DataFrame(
        y_over_h = ys_ux_trim ./ h,
        Re_u_num = real.(profile_amplitudes),
        Im_u_num = imag.(profile_amplitudes),
        Re_u_exact = real.(profile_exact),
        Im_u_exact = imag.(profile_exact)
    )

    center_exact_t = real(U_center_exact * exp.(im .* ω .* t_window))
    bulk_exact_t = real(U_bulk_exact * exp.(im .* ω .* t_window))
    shear_exact_t = real(τ_exact * exp.(im .* ω .* t_window))

    signals_df = DataFrame(
        t = t_window,
        u_center = center_window,
        u_bulk = bulk_window,
        tau_w = shear_window,
        u_center_exact = center_exact_t,
        u_bulk_exact = bulk_exact_t,
        tau_w_exact = shear_exact_t,
        force = force_window
    )

    energy_rel = abs(mean(diss_window) - mean(power_window)) / max(abs(mean(diss_window)), eps())

    return (
        profile_df = profile_df,
        signals_df = signals_df,
        summary = (
            alpha = alpha,
            omega = ω,
            dt = Δt,
            target_bulk = target_bulk,
            amp_err_center = amp_err_center,
            phase_err_center = phase_err_center,
            amp_err_bulk = amp_err_bulk,
            phase_err_bulk = phase_err_bulk,
            amp_err_shear = amp_err_shear,
            phase_err_shear = phase_err_shear,
            energy_rel = energy_rel
        ),
        profile_tag = profile_tag
    )
end

function write_outputs(results; csv_dir=nothing)
    out_dir = isnothing(csv_dir) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : csv_dir
    mkpath(out_dir)

    profiles = String[]
    signals = String[]
    for res in results
        profile_path = joinpath(out_dir, "oscillatory_channel_profile_$(res.profile_tag).csv")
        signal_path = joinpath(out_dir, "oscillatory_channel_signals_$(res.profile_tag).csv")
        CSV.write(profile_path, res.profile_df)
        CSV.write(signal_path, res.signals_df)
        push!(profiles, profile_path)
        push!(signals, signal_path)
    end

    summary_df = DataFrame(res.summary for res in results)
    summary_path = joinpath(out_dir, "oscillatory_channel_summary.csv")
    CSV.write(summary_path, summary_df)
    return (profiles=profiles, signals=signals, summary=summary_path)
end

function main(; alphas=[10.0], csv_dir=nothing)
    results = [analyze_case(α) for α in alphas]
    paths = write_outputs(results; csv_dir=csv_dir)
    return (results=results, paths=paths)
end

results = main()

@testset "Oscillatory channel Navier–Stokes" begin
    @test !isempty(results.results)
    @test all(isfile, results.paths.profiles)
    @test all(isfile, results.paths.signals)
    @test isfile(results.paths.summary)
end

println("Oscillatory channel benchmark completed. Summary: ", results.paths.summary)
