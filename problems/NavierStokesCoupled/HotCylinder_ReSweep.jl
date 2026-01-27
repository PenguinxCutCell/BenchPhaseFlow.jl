using Penguin
using LinearAlgebra
using Statistics
using Printf
using CSV
using DataFrames
using Test

"""
Sweep Reynolds numbers for a hot cylinder in channel flow.

For each inlet maximum velocity listed in `hotcylinder.txt`, solve steady
Navier-Stokes (Picard), then an unsteady advection-diffusion run with a hot
cylinder wall. Compute mean Nusselt number and compare to the reference table.
"""

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const HOTCYLINDER_DATA = joinpath(@__DIR__, "hotcylinder.txt")

function read_hotcylinder_table(path::AbstractString=HOTCYLINDER_DATA)
    rows = NamedTuple[]
    for line in readlines(path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "(") && continue
        parts = split(s)
        length(parts) < 3 && continue
        umax = parse(Float64, parts[1])
        re_d = parse(Float64, parts[2])
        nu_ref = parse(Float64, parts[3])
        push!(rows, (Umax=umax, ReD=re_d, Nu_ref=nu_ref))
    end
    return rows
end

nearest_index(vec, val) = clamp(argmin(abs.(vec .- val)), 1, length(vec))

function map_field_to_mesh(src_xs, src_ys, field, dst_xs, dst_ys)
    mapped = zeros(length(dst_xs) * length(dst_ys))
    for (j, y) in enumerate(dst_ys), (i, x) in enumerate(dst_xs)
        idx = i + (j - 1) * length(dst_xs)
        mapped[idx] = field[nearest_index(src_xs, x), nearest_index(src_ys, y)]
    end
    return mapped
end

function run_hotcylinder_case(umax;
    nx=64,
    ny=64,
    channel_length=10.0,
    channel_height=10.0,
    x0=-0.5,
    y0=-0.5,
    circle_center=(0.5, 0.0),
    circle_radius=0.2,
    mu=1.0,
    rho=1.0,
    kappa=1.0e-2,
    t_end=10.0,
    picard_tol=1e-13,
    picard_maxiter=40
)
    circle_body = (x, y, _=0.0) -> circle_radius - hypot(x - circle_center[1], y - circle_center[2])

    mesh_p = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0, y0))
    dx = mesh_p.nodes[1][2] - mesh_p.nodes[1][1]
    dy = mesh_p.nodes[2][2] - mesh_p.nodes[2][1]
    mesh_ux = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0 - 0.5 * dx, y0))
    mesh_uy = Penguin.Mesh((nx, ny), (channel_length, channel_height), (x0, y0 - 0.5 * dy))

    capacity_ux = Capacity(circle_body, mesh_ux)
    capacity_uy = Capacity(circle_body, mesh_uy)
    capacity_p = Capacity(circle_body, mesh_p)

    operator_ux = DiffusionOps(capacity_ux)
    operator_uy = DiffusionOps(capacity_uy)
    operator_p = DiffusionOps(capacity_p)

    parabolic = (x, y, t=0.0) -> begin
        xi = (y - (y0 + channel_height / 2)) / (channel_height / 2)
        umax * (1 - xi^2)
    end

    ux_left = Dirichlet((x, y, t=0.0) -> parabolic(x, y, t))
    ux_right = Dirichlet((x, y, t=0.0) -> parabolic(x, y, t))
    ux_bottom = Dirichlet((x, y, t=0.0) -> 0.0)
    ux_top = Dirichlet((x, y, t=0.0) -> 0.0)
    uy_zero = Dirichlet((x, y, t=0.0) -> 0.0)

    bc_ux = BorderConditions(Dict(
        :left => ux_left, :right => ux_right, :bottom => ux_bottom, :top => ux_top
    ))
    bc_uy = BorderConditions(Dict(
        :left => uy_zero, :right => uy_zero, :bottom => uy_zero, :top => uy_zero
    ))
    pressure_gauge = MeanPressureGauge()
    u_bc = Dirichlet(0.0)

    f_u = (x, y, z=0.0) -> 0.0
    f_p = (x, y, z=0.0) -> 0.0

    fluid = Fluid((mesh_ux, mesh_uy),
                  (capacity_ux, capacity_uy),
                  (operator_ux, operator_uy),
                  mesh_p,
                  capacity_p,
                  operator_p,
                  mu, rho, f_u, f_p)

    nu = prod(operator_ux.size)
    np = prod(operator_p.size)
    x0_vec = zeros(4 * nu + np)
    ns_solver = NavierStokesMono(fluid, (bc_ux, bc_uy), pressure_gauge, u_bc; x0=x0_vec)
    solve_NavierStokesMono_steady!(ns_solver; nlsolve_method=:picard, tol=picard_tol, maxiter=picard_maxiter, relaxation=1.0)

    uox = ns_solver.x[1:nu]
    uoy = ns_solver.x[2 * nu + 1:3 * nu]

    xs_ux, ys_ux = mesh_ux.nodes
    xs_uy, ys_uy = mesh_uy.nodes
    xp, yp = mesh_p.nodes

    Ux = reshape(uox, (length(xs_ux), length(ys_ux)))
    Uy = reshape(uoy, (length(xs_uy), length(ys_uy)))

    u_o_x = map_field_to_mesh(xs_ux, ys_ux, Ux, xp, yp)
    u_o_y = map_field_to_mesh(xs_uy, ys_uy, Uy, xp, yp)
    u_gamma = zeros(2 * length(xp) * length(yp))

    capacity_T = Capacity(circle_body, mesh_p; compute_centroids=true)
    operator_adv = ConvectionOps(capacity_T, (u_o_x, u_o_y), u_gamma)

    T_cold = 0.0
    T_hot = 1.0
    f_heat = (x, y, z, t) -> 0.0
    kappa_fn = (x, y, _=0.0) -> kappa
    phase_heat = Phase(capacity_T, operator_adv, f_heat, kappa_fn)

    ic = Dirichlet(T_hot)
    bc_b = BorderConditions(Dict(
        :left => Dirichlet(T_cold),
        :right => Dirichlet(T_cold),
        :bottom => Dirichlet(T_cold),
        :top => Dirichlet(T_cold)
    ))

    nnodes = length(xp) * length(yp)
    T0_o = fill(T_cold, nnodes)
    T0_gamma = fill(T_hot, nnodes)
    for (j, y) in enumerate(yp), (i, x) in enumerate(xp)
        if circle_body(x, y) > 0
            idx = i + (j - 1) * length(xp)
            T0_o[idx] = T_hot
        end
    end
    T0 = vcat(T0_o, T0_gamma)

    dt = 0.5 * min(dx, dy)^2
    advdiff_solver = AdvectionDiffusionUnsteadyMono(phase_heat, bc_b, ic, dt, T0, "BE")
    solve_AdvectionDiffusionUnsteadyMono!(advdiff_solver, phase_heat, dt, t_end, bc_b, ic, "CN"; method=Base.:\)

    n = prod(operator_adv.size)
    Tomega = @view advdiff_solver.states[end][1:n]
    Tgamma = @view advdiff_solver.states[end][n + 1:2 * n]

    Q = operator_adv.H' * operator_adv.Wꜝ * (operator_adv.G * Tomega + operator_adv.H * Tgamma)
    gamma_len = diag(capacity_T.Γ)

    mask_len = .!isnan.(gamma_len)
    q_n_mean = sum(Q[mask_len]) / sum(gamma_len[mask_len])
    q_n_mean *= kappa

    T_interface = sum(Tgamma[mask_len] .* gamma_len[mask_len]) / sum(gamma_len[mask_len])
    h = q_n_mean / (T_interface - T_cold)
    Nu = abs(h * (2 * circle_radius) / kappa)

    return (Nu=Nu, q_n_mean=q_n_mean, T_interface=T_interface)
end

function main(;
    nx=96,
    ny=96,
    kappa=1.0e-2,
    t_end=10.0,
    relative_tol=0.2,
    csv_path=nothing
)
    cases = read_hotcylinder_table()
    results = NamedTuple[]

    for case in cases
        res = run_hotcylinder_case(case.Umax; nx=nx, ny=ny, kappa=kappa, t_end=t_end)
        rel_err = abs(res.Nu - case.Nu_ref) / case.Nu_ref
        push!(results, (
            Umax = case.Umax,
            ReD = case.ReD,
            Nu_ref = case.Nu_ref,
            Nu = res.Nu,
            rel_err = rel_err
        ))
        @printf("Re_D = %.2f, Umax = %.2f, Nu = %.3f (ref %.3f), rel err = %.3f\n",
                case.ReD, case.Umax, res.Nu, case.Nu_ref, rel_err)
    end

    df = DataFrame(
        Umax = [r.Umax for r in results],
        ReD = [r.ReD for r in results],
        Nu_ref = [r.Nu_ref for r in results],
        Nu = [r.Nu for r in results],
        rel_err = [r.rel_err for r in results]
    )
    csv_out = isnothing(csv_path) ?
        joinpath(BENCH_ROOT, "results", "NavierStokesCoupled", "HotCylinder_ReSweep.csv") :
        csv_path
    mkpath(dirname(csv_out))
    CSV.write(csv_out, df)

    @testset "Hot cylinder Re sweep" begin
        for res in results
            @test res.rel_err ≤ relative_tol
        end
    end

    return (results = results, csv_path = csv_out, table = df)
end

results = main()
