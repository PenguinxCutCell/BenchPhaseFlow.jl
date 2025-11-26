using Penguin
using Statistics
using Printf
using CSV
using DataFrames
using Test

"""
Backward-facing step benchmark for steady incompressible Navier–Stokes (CSV only).

Solves a 2D step flow (h₂/h₁ = 2, Re = 200 default) with parabolic inlet and
no-slip walls. Reports reattachment length via wall shear, velocity profiles at
selected x/h₁ locations, bottom-wall pressure, mass-flux balance, and a small
multi-Re sweep. Outputs CSVs under `results/NavierStokes/`.
"""

###########
# Parameters
###########
const h1 = 1.0
const ratio = 2.0
const h2 = ratio * h1
const h_step = h2 - h1
const L_in = 4h1
const L_out = 20h1
const Lx = L_in + L_out
const x0 = -L_in

const ny = 80
const nx = 240

const ρ = 1.0
const U_mean = 1.0
const Re = 200.0
const ν = U_mean * h1 / Re
const μ = ρ * ν

const Um = 1.5 * U_mean  # peak for parabolic giving mean = U_mean

println("=== Backward-facing step benchmark ===")
println(@sprintf("Grid: %d × %d, Re = %.1f, ν = %.3e", nx, ny, Re, ν))

###########
# Meshes & Capacity
###########
top_y = h1
bottom_level(x) = x < 0 ? 0.0 : -h_step
step_body = (x, y, _=0.0) -> bottom_level(x) - y

mesh_p  = Penguin.Mesh((nx, ny), (Lx, top_y + h_step), (x0, -h_step))
dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
mesh_ux = Penguin.Mesh((nx, ny), (Lx, top_y + h_step), (x0 - 0.5dx, -h_step))
mesh_uy = Penguin.Mesh((nx, ny), (Lx, top_y + h_step), (x0, -h_step - 0.5dy))

capacity_ux = Capacity(step_body, mesh_ux; compute_centroids=false)
capacity_uy = Capacity(step_body, mesh_uy; compute_centroids=false)
capacity_p  = Capacity(step_body, mesh_p; compute_centroids=false)

operator_ux = DiffusionOps(capacity_ux)
operator_uy = DiffusionOps(capacity_uy)
operator_p  = DiffusionOps(capacity_p)

###########
# Boundary conditions
###########
parabolic_inlet = (x, y, t=0.0) -> begin
    if y < 0 || y > h1
        return 0.0
    end
    4Um * y * (h1 - y) / h1^2
end

no_slip = (x, y, t=0.0) -> 0.0

bc_ux = BorderConditions(Dict(
    :left   => Dirichlet(parabolic_inlet),
    :right  => Outflow(),
    :bottom => Dirichlet(no_slip),
    :top    => Dirichlet(no_slip)
))

bc_uy = BorderConditions(Dict(
    :left   => Dirichlet((x, y, t=0.0) -> 0.0),
    :right  => Outflow(),
    :bottom => Dirichlet(no_slip),
    :top    => Dirichlet(no_slip)
))

pressure_gauge = PinPressureGauge()
interface_bc = Dirichlet(0.0)

fᵤ = (x, y, z=0.0, t=0.0) -> 0.0
fₚ = (x, y, z=0.0, t=0.0) -> 0.0

fluid = Fluid((mesh_ux, mesh_uy),
              (capacity_ux, capacity_uy),
              (operator_ux, operator_uy),
              mesh_p,
              capacity_p,
              operator_p,
              μ, ρ, fᵤ, fₚ)

###########
# Solver setup
###########
nu_x = prod(operator_ux.size)
nu_y = prod(operator_uy.size)
np = prod(operator_p.size)

solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, interface_bc)

println("Solving steady Navier–Stokes...")
_, picard_iters, picard_res = solve_NavierStokesMono_steady!(
    solver;
    tol=1e-6,
    maxiter=8,
    relaxation=0.8,
    nlsolve_method=:picard,
)
println(@sprintf("Picard iterations: %d (residual %.3e)", picard_iters, picard_res))

_, newton_iters, newton_res = solve_NavierStokesMono_steady!(
    solver;
    tol=5e-9,
    maxiter=10,
    nlsolve_method=:newton,
)
println(@sprintf("Newton iterations: %d (residual %.3e)", newton_iters, newton_res))

###########
# Extract fields
###########
xs_ux = mesh_ux.nodes[1]
ys_ux = mesh_ux.nodes[2]
xs_uy = mesh_uy.nodes[1]
ys_uy = mesh_uy.nodes[2]
xs_p = mesh_p.nodes[1]
ys_p = mesh_p.nodes[2]

uωx = solver.x[1:nu_x]
uωy = solver.x[2nu_x+1:2nu_x+nu_y]
pω  = solver.x[2*(nu_x + nu_y)+1:end]

Ux = reshape(uωx, Tuple(operator_ux.size))
Uy = reshape(uωy, Tuple(operator_uy.size))
P  = reshape(pω, Tuple(operator_p.size))

trim_copy(A) = size(A, 1) > 1 && size(A, 2) > 1 ? copy(@view A[1:end-1, 1:end-1]) : copy(A)
trim_coords(v) = length(v) > 1 ? copy(v[1:end-1]) : copy(v)

Ux_trim = trim_copy(Ux)
Uy_trim = trim_copy(Uy)
P_trim  = trim_copy(P)
xs_ux_trim = trim_coords(xs_ux)
ys_ux_trim = trim_coords(ys_ux)
xs_uy_trim = trim_coords(xs_uy)
ys_uy_trim = trim_coords(ys_uy)
xs_p_trim = trim_coords(xs_p)
ys_p_trim = trim_coords(ys_p)

###########
# Diagnostics: reattachment length
###########
bottom_y(x) = x < 0 ? 0.0 : -h_step

function nearest_index(vec::AbstractVector{<:Real}, val::Real)
    clamp(argmin(abs.(vec .- val)), 1, length(vec))
end

shear_samples = Float64[]
shear_x = Float64[]
μ_inv_dy = μ / dy

for (i, x) in enumerate(xs_ux_trim)
    x < 0 && continue
    by = bottom_y(x)
    j = nearest_index(ys_ux_trim, by + dy)
    u_wall = Ux_trim[i, j]
    push!(shear_x, x)
    push!(shear_samples, μ_inv_dy * u_wall)
end

reattach_x = NaN
for k in 2:length(shear_samples)
    if shear_samples[k-1] < 0 && shear_samples[k] ≥ 0
        t = shear_samples[k-1] / (shear_samples[k-1] - shear_samples[k] + eps())
        reattach_x = shear_x[k-1] + t * (shear_x[k] - shear_x[k-1])
        break
    end
end

###########
# Velocity profiles
###########
function bilinear(xs, ys, field, xp, yp)
    x = clamp(xp, xs[1], xs[end])
    y = clamp(yp, ys[1], ys[end])
    ix_hi = searchsortedfirst(xs, x)
    ix_lo = clamp(ix_hi - 1, 1, length(xs))
    ix_hi = clamp(ix_hi, 1, length(xs))
    iy_hi = searchsortedfirst(ys, y)
    iy_lo = clamp(iy_hi - 1, 1, length(ys))
    iy_hi = clamp(iy_hi, 1, length(ys))
    x1, x2 = xs[ix_lo], xs[ix_hi]
    y1, y2 = ys[iy_lo], ys[iy_hi]
    tx = x2 ≈ x1 ? 0.0 : (x - x1) / (x2 - x1)
    ty = y2 ≈ y1 ? 0.0 : (y - y1) / (y2 - y1)
    f11 = field[ix_lo, iy_lo]
    f21 = field[ix_hi, iy_lo]
    f12 = field[ix_lo, iy_hi]
    f22 = field[ix_hi, iy_hi]
    (1 - tx) * (1 - ty) * f11 +
    tx * (1 - ty) * f21 +
    (1 - tx) * ty * f12 +
    tx * ty * f22
end

sample_positions = [2, 4, 6, 8, 10]
profile_data = Dict{Float64, Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}()

for ξ in sample_positions
    x_pos = ξ * h1
    x_pos > L_out && continue
    u_vals = Float64[]
    v_vals = Float64[]
    for y in ys_p_trim
        push!(u_vals, bilinear(xs_ux_trim, ys_ux_trim, Ux_trim, x_pos, y))
        push!(v_vals, bilinear(xs_uy_trim, ys_uy_trim, Uy_trim, x_pos, y))
    end
    profile_data[ξ] = (ys_p_trim, u_vals, v_vals)
end

###########
# Pressure along bottom wall
###########
bottom_x = Float64[]
bottom_pressure = Float64[]

for x in xs_p_trim
    x < 0 && continue
    by = bottom_y(x) + 0.5dy
    p_val = bilinear(xs_p_trim, ys_p_trim, P_trim, x, by)
    push!(bottom_x, x)
    push!(bottom_pressure, p_val)
end

###########
# Mass flux balance
###########
function integrate_face_flux(Ufield, xs, ys; column)
    values = Ufield[column, :]
    mask = (ys .>= 0) .| (xs[column] ≥ 0)
    dy_weights = diff(vcat(ys, ys[end] + dy))
    return sum(values .* dy_weights .* mask)
end

flux_in = integrate_face_flux(Ux_trim, xs_ux_trim, ys_ux_trim; column=1)
flux_out = integrate_face_flux(Ux_trim, xs_ux_trim, ys_ux_trim; column=size(Ux_trim, 1))
imbalance = flux_out - flux_in

###########
# Output CSVs
###########
out_dir = joinpath(@__DIR__, "..", "..", "results", "NavierStokes")
mkpath(out_dir)

summary_df = DataFrame(
    Re = Re,
    nu = ν,
    mu = μ,
    picard_iters = picard_iters,
    picard_residual = picard_res,
    newton_iters = newton_iters,
    newton_residual = newton_res,
    reattach_x = reattach_x,
    reattach_x_over_h1 = isnan(reattach_x) ? NaN : reattach_x / h1,
    flux_in = flux_in,
    flux_out = flux_out,
    flux_imbalance = imbalance
)
summary_path = joinpath(out_dir, "BackwardFacingStep_Summary.csv")
CSV.write(summary_path, summary_df)

profile_rows = Vector{Vector{Float64}}()
header = ["y_over_h1"]
for (ξ, _) in sort(profile_data; by=first)
    push!(header, @sprintf("u_x/h1=%.0f", ξ))
    push!(header, @sprintf("v_x/h1=%.0f", ξ))
end

max_len = maximum(length.(first.(values(profile_data))))
ys_ref = first(first(values(profile_data)))
for i in 1:max_len
    row = Float64[]
    y = i <= length(ys_ref) ? ys_ref[i] : NaN
    push!(row, y / h1)
    for (_, (y_vals, u_vals, v_vals)) in sort(profile_data; by=first)
        if i <= length(y_vals)
            push!(row, u_vals[i] / U_mean)
            push!(row, v_vals[i] / U_mean)
        else
            push!(row, NaN, NaN)
        end
    end
    push!(profile_rows, row)
end
profiles_path = joinpath(out_dir, "BackwardFacingStep_VelocityProfiles.csv")
CSV.write(profiles_path, DataFrame([header; profile_rows], :auto))

pressure_df = DataFrame(x_over_h1 = bottom_x ./ h1, p = bottom_pressure)
pressure_path = joinpath(out_dir, "BackwardFacingStep_BottomPressure.csv")
CSV.write(pressure_path, pressure_df)

shear_df = DataFrame(x_over_h1 = shear_x ./ h1, tau_w = shear_samples)
shear_path = joinpath(out_dir, "BackwardFacingStep_WallShear.csv")
CSV.write(shear_path, shear_df)

###########
# Multi-Re sweep for reattachment length (basic Picard/Newton)
###########
reference_ranges = Dict(
    100.0 => (3.5, 5.0),
    200.0 => (4.0, 6.0),
    300.0 => (5.5, 8.0),
)

function run_step_case(Re_case)
    ν_case = U_mean * h1 / Re_case
    μ_case = ρ * ν_case

    fluid_case = Fluid((mesh_ux, mesh_uy),
                       (capacity_ux, capacity_uy),
                       (operator_ux, operator_uy),
                       mesh_p,
                       capacity_p,
                       operator_p,
                       μ_case, ρ, fᵤ, fₚ)

    solver_case = NavierStokesMono(fluid_case, (bc_ux, bc_uy), pressure_gauge, interface_bc)
    solve_NavierStokesMono_steady!(solver_case; tol=1e-6, maxiter=5, relaxation=0.8, nlsolve_method=:picard)
    solve_NavierStokesMono_steady!(solver_case; tol=5e-9, maxiter=8, nlsolve_method=:newton)

    uωx_case = solver_case.x[1:nu_x]
    Ux_case = reshape(uωx_case, Tuple(operator_ux.size))
    shear = Float64[]
    xvals = Float64[]
    μ_inv = μ_case / dy
    for (i, x) in enumerate(xs_ux_trim)
        x < 0 && continue
        by = bottom_y(x)
        j = nearest_index(ys_ux_trim, by + dy)
        u_wall = Ux_case[i, j]
        isfinite(u_wall) || continue
        push!(xvals, x)
        push!(shear, μ_inv * u_wall)
    end
    reattach = NaN
    for k in 2:length(shear)
        if shear[k-1] < 0 && shear[k] ≥ 0
            t = shear[k-1] / (shear[k-1] - shear[k] + eps())
            reattach = xvals[k-1] + t * (xvals[k] - xvals[k-1])
            break
        end
    end
    return reattach
end

sweep_rows = Vector{Vector{Any}}()
for Re_case in sort(collect(keys(reference_ranges)))
    reattach = run_step_case(Re_case)
    ref_min, ref_max = reference_ranges[Re_case]
    status = isnan(reattach) ? "not-detected" :
        (ref_min <= reattach / h1 <= ref_max ? "within-range" : "outside-range")
    push!(sweep_rows, [Re_case, isnan(reattach) ? NaN : reattach / h1, ref_min, ref_max, status])
end
sweep_df = DataFrame(Re = [row[1] for row in sweep_rows],
                     Lr_over_h1 = [row[2] for row in sweep_rows],
                     ref_min = [row[3] for row in sweep_rows],
                     ref_max = [row[4] for row in sweep_rows],
                     status = [row[5] for row in sweep_rows])
sweep_path = joinpath(out_dir, "BackwardFacingStep_ReattachmentSweep.csv")
CSV.write(sweep_path, sweep_df)

println("Outputs written:")
println("  ", summary_path)
println("  ", profiles_path)
println("  ", pressure_path)
println("  ", shear_path)
println("  ", sweep_path)

###########
# Tests
###########
results_paths = (summary=summary_path, profiles=profiles_path, pressure=pressure_path, shear=shear_path, sweep=sweep_path)

@testset "Backward-facing step benchmark" begin
    @test isfile(results_paths.summary)
    @test isfile(results_paths.profiles)
    @test isfile(results_paths.pressure)
    @test isfile(results_paths.shear)
    @test isfile(results_paths.sweep)
end

println("Backward-facing step benchmark completed.")
