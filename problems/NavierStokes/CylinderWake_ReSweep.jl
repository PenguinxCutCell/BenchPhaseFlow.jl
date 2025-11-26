using Penguin
using Statistics
using Printf
using FFTW
using CSV
using DataFrames
using Test

"""
Unsteady Navier–Stokes wake benchmark around a circular cylinder (Re sweep).

Runs the cylinder wake for a small set of Reynolds numbers, extracts a probe
signal in the wake to estimate the dominant shedding frequency (Strouhal
number), and records drag/lift diagnostics. Outputs summary CSVs (no plots).
"""

###########
# Geometry and grids
###########
const nx, ny = 128, 64
const channel_length = 4.0
const channel_height = 1.0
const x0, y0 = -0.5, -0.5

const circle_center = (0.5, 0.0)
const circle_radius = 0.2
const diameter = 2 * circle_radius

circle_body = (x, y, _=0.0) -> circle_radius - hypot(x - circle_center[1], y - circle_center[2])

mesh_p  = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0, y0))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0 - 0.5 * dx, y0))
mesh_uy = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0, y0 - 0.5 * dy))

capacity_ux = Capacity(circle_body, mesh_ux)
capacity_uy = Capacity(circle_body, mesh_uy)
capacity_p  = Capacity(circle_body, mesh_p)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

###########
# Boundary conditions
###########
const Umax = 1.0
parabolic = (x, y, t=0.0) -> begin
    ξ = (y - (y0 + channel_height / 2)) / (channel_height / 2)
    Umax * (1 - ξ^2)
end

ux_left   = Dirichlet((x, y, t=0.0) -> parabolic(x, y, t))
ux_right  = Dirichlet((x, y, t=0.0) -> parabolic(x, y, t))
ux_bottom = Dirichlet((x, y, t=0.0) -> 0.0)
ux_top    = Dirichlet((x, y, t=0.0) -> 0.0)
uy_zero   = Dirichlet((x, y, t=0.0) -> 0.0)

bc_ux = BorderConditions(Dict(:left=>ux_left, :right=>Outflow(), :bottom=>Symmetry(), :top=>Symmetry()))
bc_uy = BorderConditions(Dict(:left=>uy_zero, :right=>uy_zero, :bottom=>uy_zero, :top=>uy_zero))

pressure_gauge = PinPressureGauge()
interface_bc = Dirichlet(0.0)

###########
# Helpers
###########
nearest_index(vec, val) = clamp(argmin(abs.(vec .- val)), 1, length(vec))

function estimate_strouhal_from_series(series, dt)
    s = series .- mean(series)
    N = length(s)
    if N < 8
        return (NaN, NaN)
    end
    spec = abs.(fft(s))
    freqs = (0:N-1) .* (1.0 / (dt * N))
    pos = 2:clamp(div(N, 2), 2, N)
    idx = pos[argmax(spec[pos])]
    fdom = freqs[idx]
    return fdom, spec[idx]
end

###########
# Solver factory
###########
function make_solver(μ, ρ; x0_vec=nothing)
    fluid = Fluid((mesh_ux, mesh_uy), (capacity_ux, capacity_uy), (operator_ux, operator_uy),
                  mesh_p, capacity_p, operator_p, μ, ρ, (x,y,z=0.0,t=0.0)->0.0, (x,y,z=0.0,t=0.0)->0.0)
    nu_x = prod(operator_ux.size)
    nu_y = prod(operator_uy.size)
    np = prod(operator_p.size)
    if x0_vec === nothing
        x0_vec = zeros(2 * (nu_x + nu_y) + np)
    end
    return NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, interface_bc; x0=x0_vec)
end

###########
# Run sweep
###########
function run_re_sweep(Re_list::Vector{Int})
    probe_x, probe_y = 1.0, 0.0
    xs = mesh_ux.nodes[1]
    ys = mesh_ux.nodes[2]
    ixp = nearest_index(xs, probe_x)
    iyp = nearest_index(ys, probe_y)
    LIux = LinearIndices((length(xs), length(ys)))

    nu_x = prod(operator_ux.size)
    nu_y = prod(operator_uy.size)
    np = prod(operator_p.size)

    results = Vector{NamedTuple}(undef, length(Re_list))
    series_paths = String[]

    for (idx, Re) in enumerate(Re_list)
        println("\n=== Running Re = $Re ===")
        μ = 1.0 * Umax * diameter / Re
        ρ = 1.0
        println("viscosity μ = ", μ)

        solver = make_solver(μ, ρ)

        Δt = 0.005
        T_end = Re <= 50 ? 20.0 : (Re <= 100 ? 15.0 : 10.0)

        println("Running unsteady simulation: dt=", Δt, ", T_end=", T_end)
        times, histories = solve_NavierStokesMono_unsteady!(solver; Δt=Δt, T_end=T_end, scheme=:CN)
        println("Finished: stored states = ", length(histories))

        probe_vals = Float64[]
        for hist in histories
            ux_hist = hist[1:nu_x]
            uy_hist = hist[2nu_x+1:2nu_x+nu_y]
            Uy = reshape(uy_hist, (length(xs), length(ys)))
            push!(probe_vals, Uy[ixp, iyp])
        end

        dt = times[2] - times[1]
        fdom, power = estimate_strouhal_from_series(probe_vals, dt)
        Str = fdom * diameter / Umax

        final = histories[end]
        ux_f = final[1:nu_x]
        uy_f = final[2nu_x+1:2nu_x+nu_y]
        Ux_f = reshape(ux_f, (length(xs), length(ys)))
        Uy_f = reshape(uy_f, (length(xs), length(ys)))
        speed = sqrt.(Ux_f.^2 .+ Uy_f.^2)
        max_speed = maximum(speed)
        mean_speed = mean(speed)

        ω = zeros(size(Ux_f))
        for j in 2:(size(Ux_f, 2)-1), i in 2:(size(Ux_f, 1)-1)
            dudy = (Ux_f[i, j+1] - Ux_f[i, j-1]) / (2 * dy)
            dvdx = (Uy_f[i+1, j] - Uy_f[i-1, j]) / (2 * dx)
            ω[i, j] = dvdx - dudy
        end
        max_vort = maximum(abs, ω)

        solver.x .= final
        force_diag = compute_navierstokes_force_diagnostics(solver)
        body_force = navierstokes_reaction_force_components(force_diag; acting_on=:body)
        coeffs = drag_lift_coefficients(force_diag; ρ=ρ, U_ref=Umax, length_ref=diameter, acting_on=:body)

        ref_Str = 0.2
        rel_err_Str = isfinite(Str) ? abs(Str - ref_Str) / ref_Str : NaN
        pass_strouhal = isfinite(Str) && (rel_err_Str < 0.3)
        pass_vorticity = max_vort > 0.0

        results[idx] = (
            Re = Re,
            fdom = fdom,
            Str = Str,
            rel_err_Str = rel_err_Str,
            max_speed = max_speed,
            mean_speed = mean_speed,
            max_vort = max_vort,
            Fx_body = body_force[1],
            Fy_body = body_force[2],
            Cd = coeffs.Cd,
            Cl = coeffs.Cl,
            pass_strouhal = pass_strouhal,
            pass_vorticity = pass_vorticity,
            dt = Δt,
            T_end = T_end,
            probe_power = power
        )

        # Save probe series for this Re
        series_df = DataFrame(time = times, probe_Uy = probe_vals)
        series_path = joinpath(@__DIR__, "..", "..", "results", "NavierStokes",
                               @sprintf("CylinderWake_Re%d_Probe.csv", Re))
        mkpath(dirname(series_path))
        CSV.write(series_path, series_df)
        push!(series_paths, series_path)

        println("Re = $Re summary: Str=", round(Str, sigdigits=6),
                ", rel_err_Str=", round(rel_err_Str, sigdigits=6),
                ", max_vort=", round(max_vort, sigdigits=6),
                ", Cd=", round(coeffs.Cd, sigdigits=6),
                ", Cl=", round(coeffs.Cl, sigdigits=6))
        println("checks: Str ~0.2? ", pass_strouhal, " | wake vorticity present? ", pass_vorticity)
    end

    return (results=results, series_paths=series_paths)
end

function write_summary_csv(run_results; csv_path=nothing)
    df = DataFrame(run_results.results)
    out_dir = isnothing(csv_path) ? joinpath(@__DIR__, "..", "..", "results", "NavierStokes") : dirname(csv_path)
    mkpath(out_dir)
    csv_file = isnothing(csv_path) ? joinpath(out_dir, "CylinderWake_ReSweep.csv") : csv_path
    CSV.write(csv_file, df)
    return (csv_file=csv_file, summary=df)
end

function main(; Re_list=[50, 100, 200], csv_path=nothing)
    run_results = run_re_sweep(Re_list)
    summary = write_summary_csv(run_results; csv_path=csv_path)
    return (summary=summary, series_paths=run_results.series_paths)
end

results = main()

@testset "Cylinder wake Re sweep" begin
    @test isfile(results.summary.csv_file)
    @test !isempty(results.summary.summary)
    @test all(isfile, results.series_paths)
end

println("Cylinder wake benchmark complete. Summary CSV: ", results.summary.csv_file)
