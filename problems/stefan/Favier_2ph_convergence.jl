using Penguin
using SpecialFunctions
using Roots
using CSV
using DataFrames

# Convergence study for the Favier two-phase Stefan test case.
# Outputs CSV with L2 temperature error and relative interface error + OOC.

const BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(BENCH_ROOT, "utils", "convergence.jl"))

# ------------------------------------------------------------
# Problem definition (Favier test)
# ------------------------------------------------------------
function initial_temperature(x, beta)
    return (exp(-beta * (x - 1.0)) - 1.0) / (exp(beta) - 1.0)
end

steady_state_temperature(x) = 1.0 - x
steady_state_interface = 4.0 / 5.0

function l2_volume_integrated(x, u_num, u_exact, xf)
    n = length(x)
    n == length(u_num) || error("x and u_num must have same length")
    n == length(u_exact) || error("x and u_exact must have same length")
    err = u_num .- u_exact

    total = 0.0
    for i in 1:(n - 1)
        xL = x[i]
        xR = x[i + 1]
        eL = err[i]
        eR = err[i + 1]

        if xf > xL && xf < xR
            # Split the cut cell at the interface and integrate with linear interpolation
            t = (xf - xL) / (xR - xL)
            eX = (1 - t) * eL + t * eR
            total += 0.5 * (eL^2 + eX^2) * (xf - xL)
            total += 0.5 * (eX^2 + eR^2) * (xR - xf)
        else
            total += 0.5 * (eL^2 + eR^2) * (xR - xL)
        end
    end

    volume = x[end] - x[1]
    volume > 0 || return NaN
    return sqrt(total / volume)
end

function run_favier_convergence(nx_list::Vector{Int};
    lx::Float64=1.0,
    x0::Float64=0.0,
    Tstart::Float64=0.0,
    Tend::Float64=2.0,
    alpha1::Float64=1.0,
    alpha2::Float64=1.0,
    T_bottom::Float64=1.0,
    T_top::Float64=0.0,
    T_m::Float64=1.0 / 5.0,
    rho::Float64=1.0,
    Ste::Float64=1.0,
    beta::Float64=8.041
)
    Lh = (T_bottom - T_m) / Ste

    h_vals = Float64[]
    dt_vals = Float64[]
    temp_L2 = Float64[]
    pos_rel = Float64[]

    for nx in nx_list
        lx_ext = lx + 1.0 / nx
        mesh = Penguin.Mesh((nx,), (lx_ext,), (x0,))
        x_offset = mesh.nodes[1][1] - x0
        dx = lx_ext / nx
        dt = dx^2

        xf0_phys = 1.0 / 5.0
        xf0 = xf0_phys + x_offset

        STmesh = Penguin.SpaceTimeMesh(mesh, [Tstart, Tstart + dt], tag=mesh.tag)
        body = (x, t, _=0) -> (x - xf0)
        body_c = (x, t, _=0) -> -(x - xf0)

        capacity = Capacity(body, STmesh)
        capacity_c = Capacity(body_c, STmesh)
        operator = DiffusionOps(capacity)
        operator_c = DiffusionOps(capacity_c)

        f1 = (x, y, z, t) -> 0.0
        f2 = (x, y, z, t) -> 0.0
        K1 = (x, y, z) -> alpha1
        K2 = (x, y, z) -> alpha2

        phase1 = Phase(capacity, operator, f1, K1)
        phase2 = Phase(capacity_c, operator_c, f2, K2)

        bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(
            :top => Dirichlet(T_top),
            :bottom => Dirichlet(T_bottom)
        ))
        ic = InterfaceConditions(ScalarJump(1.0, 1.0, T_m), FluxJump(1.0, 1.0, rho * Lh))

        x_phys = mesh.nodes[1] .- x_offset
        u0o1 = initial_temperature.(x_phys, beta)
        u0y1 = fill(T_m, nx + 1)
        u0o2 = initial_temperature.(x_phys, beta)
        u0y2 = fill(T_m, nx + 1)
        u0 = vcat(u0o1, u0y1, u0o2, u0y2)

        solver = MovingLiquidDiffusionUnsteadyDiph(phase1, phase2, bc_b, ic, dt, u0, mesh, "BE")
        solver, residuals, xf_log = solve_MovingLiquidDiffusionUnsteadyDiph!(
            solver, phase1, phase2, xf0, dt, Tstart, Tend, bc_b, ic, mesh, "BE";
            Newton_params=(2, 1e-10, 1e-10, 1.0),
            method=Base.:\,
            adaptive_timestep=false,
            cfl_target=0.5
        )

        xf_log_phys = xf_log .- x_offset
        xf_final = xf_log_phys[end]
        x_num = x_phys

        u1_bulk = solver.x[1:(nx + 1)]
        u2_bulk = solver.x[2 * (nx + 1) + 1:3 * (nx + 1)]
        u_num = similar(x_num)
        mask_left = x_num .<= xf_final
        u_num[mask_left] .= u1_bulk[mask_left]
        u_num[.!mask_left] .= u2_bulk[.!mask_left]

        u_steady = steady_state_temperature.(x_num)
        L2_error = l2_volume_integrated(x_num, u_num, u_steady, xf_final)

        pos_err = abs(xf_final - steady_state_interface)
        rel_pos_err = pos_err / (abs(steady_state_interface) > 0 ? abs(steady_state_interface) : eps())

        push!(h_vals, dx)
        push!(dt_vals, dt)
        push!(temp_L2, L2_error)
        push!(pos_rel, rel_pos_err)
    end

    ooc_temp = compute_pairwise_orders(h_vals, temp_L2)
    ooc_pos = compute_pairwise_orders(h_vals, pos_rel)

    return (
        h_vals=h_vals,
        dt_vals=dt_vals,
        temp_L2=temp_L2,
        pos_rel=pos_rel,
        ooc_temp=ooc_temp,
        ooc_pos=ooc_pos
    )
end

function write_favier_csv(nx_list, data; csv_path=nothing)
    results_dir = isnothing(csv_path) ? joinpath(BENCH_ROOT, "results", "stefan") : dirname(csv_path)
    mkpath(results_dir)
    csv_out = isnothing(csv_path) ? joinpath(results_dir, "Favier_2ph_convergence.csv") : csv_path

    df = DataFrame(
        nx = nx_list,
        h = data.h_vals,
        dt = data.dt_vals,
        L2_error = data.temp_L2,
        rel_pos_error = data.pos_rel,
        ooc_L2 = data.ooc_temp,
        ooc_pos = data.ooc_pos
    )
    CSV.write(csv_out, df)
    return (csv_path=csv_out, table=df)
end

function main(; nx_list=nothing, csv_path=nothing)
    nx_vals = isnothing(nx_list) ? [4, 8, 16, 32, 64, 128] : nx_list
    data = run_favier_convergence(nx_vals)
    csv_info = write_favier_csv(nx_vals, data; csv_path=csv_path)
    return (data=data, csv_path=csv_info.csv_path, table=csv_info.table)
end

results = main()
println("Favier 2ph convergence written to:\n  ", results.csv_path)
