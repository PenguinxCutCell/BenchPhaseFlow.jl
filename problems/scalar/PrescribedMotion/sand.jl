using CairoMakie

"""
Exact-solution snapshots for the 2-phase vertically oscillating interface
benchmark (Heat_2ph_2D_MovingVertical).
"""

struct MovingVerticalParams
    lx::Float64
    ly::Float64
    x0::Float64
    y0::Float64
    Tend::Float64
    D_plus::Float64
    D_minus::Float64
    cp_plus::Float64
    cp_minus::Float64
    s0::Float64
    A::Float64
    ω::Float64
end

MovingVerticalParams(; lx=4.0, ly=4.0, x0=0.0, y0=0.0, Tend=0.1,
                    D_plus=1.0, D_minus=1.0, cp_plus=1.0, cp_minus=1.0,
                    s0=2.0, A=1.0, ω=4π) = MovingVerticalParams(
    lx, ly, x0, y0, Tend, D_plus, D_minus, cp_plus, cp_minus, s0, A, ω
)

s_val(t, p::MovingVerticalParams) = p.s0 + p.A * sin(p.ω * t)

g1(x, p::MovingVerticalParams) = (x - p.x0) * (p.lx - (x - p.x0))
g2(y, p::MovingVerticalParams) = (y - p.y0) * (p.ly - (y - p.y0))
g(x, y, p::MovingVerticalParams) = g1(x, p) * g2(y, p)
θ(t) = exp(-t)

function u1_exact(p::MovingVerticalParams)
    return (x, y, t) -> x >= s_val(t, p) ? p.D_minus * (x - s_val(t, p)) * g(x, y, p) * θ(t) : 0.0
end

function u2_exact(p::MovingVerticalParams)
    return (x, y, t) -> x < s_val(t, p) ? p.D_plus * (x - s_val(t, p)) * g(x, y, p) * θ(t) : 0.0
end

# Assemble the piecewise exact solution over a grid.
function snapshot_exact(p::MovingVerticalParams, t; nx=200, ny=200)
    xs = range(p.x0; stop=p.x0 + p.lx, length=nx)
    ys = range(p.y0; stop=p.y0 + p.ly, length=ny)
    u1 = u1_exact(p)
    u2 = u2_exact(p)
    field = Matrix{Float64}(undef, ny, nx)
    s = s_val(t, p)
    for (j, y) in enumerate(ys)
        for (i, x) in enumerate(xs)
            field[j, i] = x >= s ? u1(x, y, t) : u2(x, y, t)
        end
    end
    return xs, ys, field, s
end

function plot_exact_snapshots(; ts=(0.0, 0.125, 0.25, 0.375),
    params::MovingVerticalParams=MovingVerticalParams(),
    output_path=joinpath(@__DIR__, "moving_vertical_exact_snapshots.png"),
    nx=200, ny=200)

    snapshots = [snapshot_exact(params, t; nx=nx, ny=ny) for t in ts]
    maxabs = maximum(abs, vcat([vec(f) for (_, _, f, _) in snapshots]...))
    colorrange = (-maxabs, maxabs)

    fig = Figure(size=(1400, 900))
    color_ref = nothing

    for (idx, (t, (xs, ys, field, s))) in enumerate(zip(ts, snapshots))
        row = cld(idx, 2)
        col = idx - 2 * (row - 1)
        ax = Axis(fig[row, col]; title="t = $(round(t, digits=4))", aspect=DataAspect(), xlabel="x", ylabel="y")
        xs_vec = collect(xs)
        ys_vec = collect(ys)
        
        hm = heatmap!(ax, xs_vec, ys_vec, field'; colormap=:balance, colorrange=colorrange)
        lines!(ax, [s, s], [minimum(ys_vec), maximum(ys_vec)];
               color=:red, linestyle=:dash, linewidth=2,
               label=idx == 1 ? "interface" : "")
        color_ref === nothing && (color_ref = hm)
        xlims!(ax, first(xs_vec), last(xs_vec))
        ylims!(ax, first(ys_vec), last(ys_vec))
        idx == 1 && axislegend(ax; position=:lt)
    end

    Colorbar(fig[:, end + 1], color_ref; label="Numerical value")
    save(output_path, fig)
    return output_path
end

function main()
    out = plot_exact_snapshots()
    println("Saved exact-solution snapshots to: $(out)")
end

main()
