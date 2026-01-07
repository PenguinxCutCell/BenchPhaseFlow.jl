using LsqFit
using DataFrames
using Penguin

"""
    count_inside_cells(capacity)

Return the number of fully inside cells (cell type == 1).
"""
count_inside_cells(capacity) = count(x -> x == 1, capacity.cell_types)

# Keep only indices that are at least `trim_layers` cells away from the boundary.
# Returns `nothing` when trimming is disabled or not feasible.
function build_interior_mask(capacity::Capacity, trim_layers::Int)
    trim_layers <= 0 && return nothing

    mesh = capacity.mesh
    dims = ntuple(i -> length(mesh.nodes[i]), length(mesh.nodes))

    if any(d -> d <= 2 * trim_layers, dims)
        @warn "Cannot trim $trim_layers layers: mesh too small (dims = $dims)"
        return nothing
    end

    mask = falses(prod(dims))
    li = LinearIndices(dims)
    ranges = ntuple(i -> (1 + trim_layers):(dims[i] - trim_layers), length(dims))
    for ci in CartesianIndices(ranges)
        mask[li[ci]] = true
    end

    if length(mask) != length(capacity.cell_types)
        @warn "Trim mask length $(length(mask)) does not match cell_types length $(length(capacity.cell_types)); skipping trim"
        return nothing
    end

    return mask
end

apply_trim(idx::Vector{Int}, mask) = mask === nothing ? idx : [i for i in idx if mask[i]]

"""
    check_h1_convergence(grad_analytical, solver, capacity::Capacity{D}, operator; p::Real=2, relative::Bool=false) where D

Compute gradient-based error norms (H1 semi-norm) between numerical and analytical
solutions. `grad_analytical` can be either a tuple of `D` component functions or a
single function returning `D` gradient components.

Returns analytical gradients, numerical gradients, and the norms over all, full,
cut, and empty cells (mirroring `check_convergence`).
"""
function check_h1_convergence(grad_analytical, solver, capacity::Capacity{D}, operator; p::Real=2, relative::Bool=false, trim_layers::Int=2) where D
    cell_centroids = capacity.C_ω
    n_cells = length(cell_centroids)

    grad_num_flat = ∇(operator, solver.x)
    expected_len = D * n_cells
    expected_len == length(grad_num_flat) || error("Expected gradient of length $expected_len, got $(length(grad_num_flat))")

    comp_len = n_cells
    grad_num = [grad_num_flat[(i-1)*comp_len + 1 : i*comp_len] for i in 1:D]
    grad_ana = [zeros(Float64, n_cells) for _ in 1:D]

    for (idx, c) in enumerate(cell_centroids)
        coords = ntuple(i -> c[i], D)
        raw_vals = grad_analytical isa NTuple{D,Function} ?
            ntuple(i -> grad_analytical[i](coords...), D) :
            grad_analytical(coords...)
        vals = Tuple(raw_vals)
        length(vals) == D || error("Analytical gradient must provide $D components, got $(length(vals))")
        for d in 1:D
            grad_ana[d][idx] = vals[d]
        end
    end

    grad_err = [grad_ana[d] .- grad_num[d] for d in 1:D]
    err_mag = [sqrt(sum(abs2(grad_err[d][i]) for d in 1:D)) for i in 1:n_cells]
    ana_mag = [sqrt(sum(abs2(grad_ana[d][i]) for d in 1:D)) for i in 1:n_cells]

    cell_types = capacity.cell_types
    mask = build_interior_mask(capacity, trim_layers)
    idx_all = apply_trim(findall((cell_types .== 1) .| (cell_types .== -1)), mask)
    idx_full = apply_trim(findall(cell_types .== 1), mask)
    idx_cut = apply_trim(findall(cell_types .== -1), mask)
    idx_empty = apply_trim(findall(cell_types .== 0), mask)

    if relative
        global_err = relative_lp_norm(err_mag, idx_all, p, capacity, ana_mag)
        full_err = relative_lp_norm(err_mag, idx_full, p, capacity, ana_mag)
        cut_err = relative_lp_norm(err_mag, idx_cut, p, capacity, ana_mag)
        empty_err = relative_lp_norm(err_mag, idx_empty, p, capacity, ana_mag)
    else
        global_err = Penguin.lp_norm(err_mag, idx_all, p, capacity)
        full_err = Penguin.lp_norm(err_mag, idx_full, p, capacity)
        cut_err = Penguin.lp_norm(err_mag, idx_cut, p, capacity)
        empty_err = Penguin.lp_norm(err_mag, idx_empty, p, capacity)
    end

    println("All cells H1 L$p norm   = $global_err")
    println("Full cells H1 L$p norm  = $full_err")
    println("Cut cells H1 L$p norm   = $cut_err")
    println("Empty cells H1 L$p norm = $empty_err")

    return (grad_ana, grad_num, global_err, full_err, cut_err, empty_err)
end

"""
    check_convergence_diph(u1_analytical, u2_analytical, solver,
                           capacity1::Capacity{D}, capacity2::Capacity{D};
                           p::Real=2, relative::Bool=false, trim_layers::Int=2) where D

Local override of Penguin's diphasic convergence helper with dimension-agnostic
indexing. Uses the length of each capacity rather than 1D-specific assumptions.
Falls back to `solver.x` when no state history is stored.
"""
function check_convergence_diphh(u1_analytical::Function, u2_analytical::Function, solver,
                                capacity1::Capacity{D}, capacity2::Capacity{D};
                                p::Real=2, relative::Bool=false, trim_layers::Int=2) where D
    # Extract centroids
    c1 = capacity1.C_ω
    c2 = capacity2.C_ω
    n1 = length(c1)
    n2 = length(c2)
    n1 == n2 || @warn "Capacity sizes differ: n1=$n1, n2=$n2; using separate lengths for extraction."

    # Evaluate analytical fields
    u1_ana = if D == 1
        map(c -> u1_analytical(c[1]), c1)
    elseif D == 2
        map(c -> u1_analytical(c[1], c[2]), c1)
    elseif D == 3
        map(c -> u1_analytical(c[1], c[2], c[3]), c1)
    elseif D == 4
        map(c -> u1_analytical(c[1], c[2], c[3], c[4]), c1)
    else
        error("Unsupported dimension: $D")
    end

    u2_ana = if D == 1
        map(c -> u2_analytical(c[1]), c2)
    elseif D == 2
        map(c -> u2_analytical(c[1], c[2]), c2)
    elseif D == 3
        map(c -> u2_analytical(c[1], c[2], c[3]), c2)
    elseif D == 4
        map(c -> u2_analytical(c[1], c[2], c[3], c[4]), c2)
    else
        error("Unsupported dimension: $D")
    end

    x = !isempty(solver.states) ? solver.states[end] : solver.x
    x === nothing && error("Solver state unavailable: solver.states is empty and solver.x is nothing")

    expected = 2n1 + 2n2
    length(x) >= expected || error("Solver state too small for diphasic extraction: length(x)=$(length(x)) < $expected")

    # Layout: [uω1, uγ1, uω2, uγ2]
    u1_num = x[1:n1]
    u2_num = x[2n1 + 1 : 2n1 + n2]

    err1 = u1_ana .- u1_num
    err2 = u2_ana .- u2_num

    ctype1 = capacity1.cell_types
    ctype2 = capacity2.cell_types

    mask1 = build_interior_mask(capacity1, trim_layers)
    mask2 = build_interior_mask(capacity2, trim_layers)

    idx_all1 = apply_trim(findall((ctype1 .== 1) .| (ctype1 .== -1)), mask1)
    idx_full1 = apply_trim(findall(ctype1 .== 1), mask1)
    idx_cut1 = apply_trim(findall(ctype1 .== -1), mask1)
    idx_empty1 = apply_trim(findall(ctype1 .== 0), mask1)

    idx_all2 = apply_trim(findall((ctype2 .== 1) .| (ctype2 .== -1)), mask2)
    idx_full2 = apply_trim(findall(ctype2 .== 1), mask2)
    idx_cut2 = apply_trim(findall(ctype2 .== -1), mask2)
    idx_empty2 = apply_trim(findall(ctype2 .== 0), mask2)

    # Phase 1 norms
    if relative
        global_err1 = Penguin.relative_lp_norm(err1, idx_all1, p, capacity1, u1_ana)
        full_err1 = Penguin.relative_lp_norm(err1, idx_full1, p, capacity1, u1_ana)
        cut_err1 = Penguin.relative_lp_norm(err1, idx_cut1, p, capacity1, u1_ana)
        empty_err1 = Penguin.relative_lp_norm(err1, idx_empty1, p, capacity1, u1_ana)
    else
        global_err1 = Penguin.lp_norm(err1, idx_all1, p, capacity1)
        full_err1 = Penguin.lp_norm(err1, idx_full1, p, capacity1)
        cut_err1 = Penguin.lp_norm(err1, idx_cut1, p, capacity1)
        empty_err1 = Penguin.lp_norm(err1, idx_empty1, p, capacity1)
    end

    # Phase 2 norms
    if relative
        global_err2 = Penguin.relative_lp_norm(err2, idx_all2, p, capacity2, u2_ana)
        full_err2 = Penguin.relative_lp_norm(err2, idx_full2, p, capacity2, u2_ana)
        cut_err2 = Penguin.relative_lp_norm(err2, idx_cut2, p, capacity2, u2_ana)
        empty_err2 = Penguin.relative_lp_norm(err2, idx_empty2, p, capacity2, u2_ana)
    else
        global_err2 = Penguin.lp_norm(err2, idx_all2, p, capacity2)
        full_err2 = Penguin.lp_norm(err2, idx_full2, p, capacity2)
        cut_err2 = Penguin.lp_norm(err2, idx_cut2, p, capacity2)
        empty_err2 = Penguin.lp_norm(err2, idx_empty2, p, capacity2)
    end

    global_err = max(global_err1, global_err2)
    full_err = max(full_err1, full_err2)
    cut_err = max(cut_err1, cut_err2)
    empty_err = max(empty_err1, empty_err2)

    println("\n=== Phase 1 Errors ===")
    println("Phase 1 - All cells L$p norm   = $global_err1")
    println("Phase 1 - Full cells L$p norm  = $full_err1")
    println("Phase 1 - Cut cells L$p norm   = $cut_err1")
    println("Phase 1 - Empty cells L$p norm = $empty_err1")

    println("\n=== Phase 2 Errors ===")
    println("Phase 2 - All cells L$p norm   = $global_err2")
    println("Phase 2 - Full cells L$p norm  = $full_err2")
    println("Phase 2 - Cut cells L$p norm   = $cut_err2")
    println("Phase 2 - Empty cells L$p norm = $empty_err2")

    println("\n=== Combined Errors (maximum of both phases) ===")
    println("Combined - All cells L$p norm   = $global_err")
    println("Combined - Full cells L$p norm  = $full_err")
    println("Combined - Cut cells L$p norm   = $cut_err")
    println("Combined - Empty cells L$p norm = $empty_err")

    return (
        (u1_ana, u2_ana),
        (u1_num, u2_num),
        (global_err1, global_err2, global_err),
        (full_err1, full_err2, full_err),
        (cut_err1, cut_err2, cut_err),
        (empty_err1, empty_err2, empty_err)
    )
end

# Backward-compatibility alias for older scripts that used the misspelled name.
check_convergence_diphh(args...; kwargs...) = check_convergence_diph(args...; kwargs...)


"""
    compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals; use_last=3)

Return rounded convergence-rate estimates for the global, full-cell, and cut-cell
errors. Fits are performed on `log(err)` vs `log(h)` using Linear Least Squares.
"""
function compute_orders(h_vals, err_vals, err_full_vals, err_cut_vals)
    function fit_model(x, p)
        p[1] .* x .+ p[2]
    end

    function safe_fit(h, err, use_last_n)
        mask = err .> 0
        if count(mask) < 2
            return NaN
        end
        h_pos = h[mask]
        err_pos = err[mask]
        log_h = log.(h_pos)
        log_err = log.(err_pos)

        n = min(use_last_n, length(log_h))
        idx = length(log_h) - n + 1 : length(log_h)
        fit_result = curve_fit(fit_model, log_h[idx], log_err[idx], [-1.0, 0.0])
        return fit_result.param[1]
    end

    p_all_all  = safe_fit(h_vals, err_vals, length(h_vals))
    p_full_all = safe_fit(h_vals, err_full_vals, length(h_vals))
    p_cut_all  = safe_fit(h_vals, err_cut_vals, length(h_vals))

    p_all  = safe_fit(h_vals, err_vals, 3)
    p_full = safe_fit(h_vals, err_full_vals, 3)
    p_cut  = safe_fit(h_vals, err_cut_vals, 3)

    round_or_nan(x) = isnan(x) ? NaN : round(x, digits=1)

    return (
        all = round_or_nan(p_all),
        full = round_or_nan(p_full),
        cut = round_or_nan(p_cut),
        all_all = round_or_nan(p_all_all),
        full_all = round_or_nan(p_full_all),
        cut_all = round_or_nan(p_cut_all)
    )
end

"""
    compute_pairwise_orders(h_vals, err_vals)

Return per-level convergence orders computed between consecutive mesh levels:
`p_i = log(e_{i-1}/e_i) / log(h_{i-1}/h_i)` for `i >= 2`. Entries with
insufficient data are filled with `NaN`.
"""
function compute_pairwise_orders(h_vals, err_vals)
    n = length(h_vals)
    orders = fill(NaN, n)
    for i in 2:n
        h_prev = h_vals[i-1]
        h_curr = h_vals[i]
        e_prev = err_vals[i-1]
        e_curr = err_vals[i]
        if h_prev > 0 && h_curr > 0 && e_prev > 0 && e_curr > 0 && h_prev != h_curr
            orders[i] = log(e_prev / e_curr) / log(h_prev / h_curr)
        end
    end
    return orders
end

"""
    make_convergence_dataframe(method_name, data)

Create a `DataFrame` with convergence information for a given method.
"""
function make_convergence_dataframe(method_name, data)
    n = length(data.h_vals)
    lp_label = Vector{Union{Missing,String}}(undef, n)
    norm_value = haskey(data, :norm) ? data.norm : nothing

    if isnothing(norm_value)
        fill!(lp_label, missing)
    else
        fill!(lp_label, "L^$(norm_value)")
    end

    pair_all = compute_pairwise_orders(data.h_vals, data.err_vals)
    pair_full = compute_pairwise_orders(data.h_vals, data.err_full_vals)
    pair_cut = compute_pairwise_orders(data.h_vals, data.err_cut_vals)

    df = DataFrame(
        method = fill(method_name, n),
        h = data.h_vals,
        lp_norm = lp_label,
        inside_cells = data.inside_cells,
        inside_cells_by_dim = data.inside_cells_by_dim,
        all_err = data.err_vals,
        full_err = data.err_full_vals,
        cut_err = data.err_cut_vals,
        empty_err = data.err_empty_vals,
        pair_order_all = pair_all,
        pair_order_full = pair_full,
        pair_order_cut = pair_cut
    )

    return df
end

"""
    make_diphasic_convergence_dataframe(method_name, data)

Create a `DataFrame` with per-phase and combined convergence information for a
two-phase configuration. `data` must provide the same fields as the monophasic
case plus per-phase errors and cell counts.
"""
function make_diphasic_convergence_dataframe(method_name, data)
    n = length(data.h_vals)
    lp_label = Vector{Union{Missing,String}}(undef, n)
    norm_value = haskey(data, :norm) ? data.norm : nothing

    if isnothing(norm_value)
        fill!(lp_label, missing)
    else
        fill!(lp_label, "L^$(norm_value)")
    end

    pair_all = compute_pairwise_orders(data.h_vals, data.err_vals)
    pair_full = compute_pairwise_orders(data.h_vals, data.err_full_vals)
    pair_cut = compute_pairwise_orders(data.h_vals, data.err_cut_vals)

    df = DataFrame(
        method = fill(method_name, n),
        h = data.h_vals,
        lp_norm = lp_label,
        inside_cells = data.inside_cells,
        inside_cells_phase1 = data.inside_cells_phase1,
        inside_cells_phase2 = data.inside_cells_phase2,
        inside_cells_by_dim = data.inside_cells_by_dim,
        all_err = data.err_vals,
        full_err = data.err_full_vals,
        cut_err = data.err_cut_vals,
        empty_err = data.err_empty_vals,
        phase1_all_err = data.phase1_all_errs,
        phase1_full_err = data.phase1_full_errs,
        phase1_cut_err = data.phase1_cut_errs,
        phase1_empty_err = data.phase1_empty_errs,
        phase2_all_err = data.phase2_all_errs,
        phase2_full_err = data.phase2_full_errs,
        phase2_cut_err = data.phase2_cut_errs,
        phase2_empty_err = data.phase2_empty_errs,
        pair_order_all = pair_all,
        pair_order_full = pair_full,
        pair_order_cut = pair_cut
    )

    return df
end
