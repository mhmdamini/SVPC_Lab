# =============================================================================
# Reformulations used by the cut generator.
#
#   standard form :  max  c'x          s.t.  A x <= b,         x >= 0
#   normal form   :  max  c'x + 0's    s.t.  A x + I s == b,   x, s >= 0
#
# The normal form adds one slack per row so that the simplex tableau (and
# therefore the cone rays at the optimal vertex) can be read off directly.
# =============================================================================

"""
    create_standard_form_model(model) -> Model

Rewrite an LP as `max c'x  s.t.  A x <= b, x >= 0`:

* `>=` rows are negated, equality rows become a `<=`/`>=` pair;
* finite variable upper bounds and positive lower bounds become explicit rows;
* a `Min` objective is turned into `Max` of its negation.

Variables keep their order (and names), so a solution of the standard form
maps one-to-one onto the original variables. All variables are assumed to be
non-negative. Each row is named after its origin (`"c1"`, `"c1 (≥)"` for the
negated half of a range/equality row, `"x1 ≤ 1"` for a bound) so tableaux can
be labelled.
"""
function create_standard_form_model(model::Model)
    data = lp_matrix_data(model)
    A = Matrix(data.A)
    n = length(data.variables)
    var_names = [isempty(name(v)) ? "x$j" : name(v) for (j, v) in enumerate(data.variables)]

    std_model = Model()
    @variable(std_model, x[1:n] >= 0)
    set_name.(x, var_names)
    add_row!(expr, rhs, label) = set_name(@constraint(std_model, expr <= rhs), label)

    # Linear constraints  b_lower <= A x <= b_upper  →  "<=" rows only.
    for i in 1:size(A, 1)
        lo, up = data.b_lower[i], data.b_upper[i]
        label = name(data.affine_constraints[i])
        label = isempty(label) ? "row $i" : label
        both = isfinite(lo) && isfinite(up)
        isfinite(up) && add_row!(A[i, :]' * x, up, both ? "$label (≤)" : label)
        isfinite(lo) && add_row!(-A[i, :]' * x, -lo, both ? "$label (≥)" : label)
    end

    # Variable bounds that are not already implied by x >= 0.
    for j in 1:n
        lo, up = data.x_lower[j], data.x_upper[j]
        isfinite(up) && add_row!(1.0 * x[j], up, "$(var_names[j]) ≤ $(number_text(up))")
        isfinite(lo) && lo > 0 && add_row!(-1.0 * x[j], -lo, "$(var_names[j]) ≥ $(number_text(lo))")
    end

    # Always maximise; a minimisation objective is negated.
    sign = data.sense == MOI.MIN_SENSE ? -1 : 1
    @objective(std_model, Max, sign * (data.c' * x + data.c_offset))
    return std_model
end

"""Short text for a number in a label: `1` rather than `1.0`."""
number_text(x::Real) = isinteger(x) ? string(Int(x)) : string(x)

"""
    create_normal_form_model(model) -> Model

Add one slack variable per row of a standard-form model, giving
`max c'x  s.t.  [A I] [x; s] == b, x, s >= 0`. The structural variables are
named `β[i]`, the slacks `s[i]`, and constraint `i` is named `"i"` so that the
basis status of each row can be mapped back to its slack later.
"""
function create_normal_form_model(model::Model)
    data = lp_matrix_data(model)
    A = Matrix(data.A)
    m, n = size(A)
    c = vcat(data.c, zeros(m))

    norm_model = Model()
    @variable(norm_model, β[1:n] >= 0)
    @variable(norm_model, s[1:m] >= 0)
    xs = vcat(β, s)
    rows = @constraint(norm_model, hcat(A, I(m)) * xs .== data.b_upper)
    @objective(norm_model, Max, c' * xs + data.c_offset)

    for (i, con) in enumerate(rows)
        set_name(con, "$i")
    end
    return norm_model
end
