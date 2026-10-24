# =============================================================================
# V-polyhedral disjunctive cuts (VPCs).
#
# For a split disjunction  x_k <= ⌊x̄_k⌋  ∨  x_k >= ⌈x̄_k⌉  each side t gives an
# LP whose optimal vertex p_t and cone rays R_t describe a relaxation of that
# side. Any (α, β) with
#
#     α'p_t >= β     and     α'r >= 0   for every ray r ∈ R_t,  t = 1, 2
#
# is valid for the union of both sides, hence for every integer point. The
# point-ray LP (PRLP) picks the deepest such cut with respect to the current
# LP optimum ω, using the normalisation β = 1.
# =============================================================================

"""
    find_fractional_elements(values; atol = 1e-6) -> (fractional_values, indices)

Entries of `values` that are farther than `atol` from the nearest integer,
together with their positions.
"""
function find_fractional_elements(values::AbstractVector{<:Real}; atol::Real = 1e-6)
    indices = findall(v -> abs(v - round(v)) > atol, values)
    return Float64.(values[indices]), indices
end

"""
    fractional_integer_variables(relaxed_model, integer_vars; atol = 1e-6) -> Vector{Int}

Indices of the integer variables that take a fractional value in the optimal
solution of the (solved) LP relaxation. Empty means the LP optimum is integral.
"""
function fractional_integer_variables(relaxed_model::Model, integer_vars; atol::Real = 1e-6)
    values, _ = get_variable_values(relaxed_model)
    _, fractional = find_fractional_elements(values; atol)
    return intersect(fractional, integer_vars)
end

"""
    create_disjunctive_subproblems(model, var_index, integer_vars;
                                   names = ("Disjunction1", "Disjunction2"),
                                   optimizer = HiGHS.Optimizer)
        -> (subproblems, D0)

Solve `model` (in standard form) and split on variable `var_index`, which must
be an integer variable with a fractional value x̄. Returns two copies of the
model with the extra rows

    side 1:   x_k <= ⌊x̄⌋
    side 2:  -x_k <= -(⌊x̄⌋ + 1)

and `D0 = [⌊x̄⌋, -(⌊x̄⌋ + 1)]`, the right-hand sides of these rows. The split
row is the last row of each side (`monoidal_strengthening` relies on this).
"""
function create_disjunctive_subproblems(model::Model, var_index::Int, integer_vars;
                                        names = ("Disjunction1", "Disjunction2"),
                                        optimizer = HiGHS.Optimizer)
    var_index in integer_vars || error("Variable $var_index is not an integer variable.")

    optimize_model!(model; optimizer)
    values, _ = get_variable_values(model)
    x̄ = values[var_index]
    abs(x̄ - round(x̄)) > 1e-6 || error("Variable $var_index has integral value $x̄; nothing to split on.")
    floor_x̄ = floor(x̄)

    side1 = copy(model)
    x1 = all_variables(side1)[var_index]
    set_name(@constraint(side1, x1 <= floor_x̄), names[1])

    side2 = copy(model)
    x2 = all_variables(side2)[var_index]
    set_name(@constraint(side2, -x2 <= -(floor_x̄ + 1)), names[2])

    return [side1, side2], [floor_x̄, -(floor_x̄ + 1)]
end

"""
    process_subproblem(subproblem; optimizer = HiGHS.Optimizer)

Convert one side of the disjunction to normal form, solve it and read off its
simplex information. Returns a named tuple with

* `first_tableau`, `optimal_tableau` – initial and optimal tableaux (`DataFrame`s);
* `rays`        – cone rays at the optimal vertex (one column per ray);
* `point`       – optimal vertex, structural variables only;
* `nonbasic`    – names of the non-basic variables (one per ray);
* `normal_model` – the solved normal-form model.
"""
function process_subproblem(subproblem::Model; optimizer = HiGHS.Optimizer)
    normal_model = create_normal_form_model(subproblem)
    optimize_model!(normal_model; optimizer)

    first_tableau = generate_first_simplex_tableau(normal_model)
    opt = generate_optimal_simplex_tableau(normal_model)
    rays = compute_rays(opt.tableau, opt.structural, opt.nonbasic)
    values, _ = get_variable_values(normal_model)

    return (first_tableau = first_tableau, optimal_tableau = opt.tableau, rays = rays,
            point = values[1:length(opt.structural)], nonbasic = opt.nonbasic,
            normal_model = normal_model)
end

"""
    solve_prlp(ω, points, rays; normalization = :beta, optimizer = HiGHS.Optimizer)
        -> (prlp, α, β)

Point-ray LP: find a cut `α'x >= β` valid for every point and ray,

    α'p_t >= β        for every point p_t
    α'r   >= 0        for every ray r (columns of rays[t])

that is as violated as possible at `ω`, the LP optimum to be cut off.
`normalization` selects how the cut is scaled:

* `:beta` – fix `β = 1` and minimise `α'ω` (the original formulation). Unbounded
  when a valid cut with `β = 0` separates `ω`, infeasible when every valid
  cut has `β <= 0`.
* `:l1`   – keep `β` free, require `‖α‖₁ <= 1` and minimise `α'ω - β`. Always
  bounded; returns the most violated cut under this scaling.
* `:auto` – try `:beta`, fall back to `:l1` if that LP has no optimal solution
  (announced on stdout unless `verbose = false`). The normalisation actually
  used can be read back with `is_fixed(prlp[:β])` (`true` ⇔ `:beta`).
"""
function solve_prlp(ω::AbstractVector, points, rays; normalization::Symbol = :beta,
                    optimizer = HiGHS.Optimizer, verbose::Bool = true)
    if normalization == :auto
        try
            return solve_prlp(ω, points, rays; normalization = :beta, optimizer)
        catch err
            err isa ErrorException && startswith(err.msg, "PRLP did not solve") || rethrow()
            verbose && println("PRLP with β = 1 has no optimal solution; using the ℓ1-normalised PRLP instead.")
            return solve_prlp(ω, points, rays; normalization = :l1, optimizer)
        end
    end
    normalization in (:beta, :l1) || error("Unknown normalization :$normalization (use :beta, :l1 or :auto).")

    n = length(ω)
    prlp = Model()
    @variable(prlp, α[1:n])
    @variable(prlp, β)
    for p in points
        @constraint(prlp, p[1:n]' * α >= β)
    end
    for R in rays, j in axes(R, 2)
        @constraint(prlp, R[:, j]' * α >= 0)
    end

    if normalization == :beta
        fix(β, 1.0)
        @objective(prlp, Min, ω' * α)
    else
        @variable(prlp, abs_α[1:n] >= 0)          # abs_α ≥ |α| componentwise
        @constraint(prlp, abs_α .>= α)
        @constraint(prlp, abs_α .>= -α)
        @constraint(prlp, sum(abs_α) <= 1)
        @objective(prlp, Min, ω' * α - β)
    end

    optimize_model!(prlp; optimizer)
    status = termination_status(prlp)
    status == MOI.OPTIMAL || error("PRLP did not solve to optimality (status: $status).")
    return prlp, value.(α), value(β)
end

"""
    add_cut(model, α, β; verbose = true) -> Model

Return a copy of `model` with the cut `α'x >= β` added (`x = all_variables`).
"""
function add_cut(model::Model, α::AbstractVector, β::Real; verbose::Bool = true)
    cut_model = copy(model)
    cut = @constraint(cut_model, α' * all_variables(cut_model) >= β)
    set_name(cut, "cut")
    verbose && println("Added cut: ", cut)
    return cut_model
end
