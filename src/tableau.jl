# =============================================================================
# Simplex tableaux of a solved normal-form model, and the extreme rays of the
# LP cone at its optimal vertex.
#
# Tableau convention (one row per basic variable plus an objective row "z"):
#
#     x_B + B⁻¹N x_N = B⁻¹b
#     z  + (c_B'B⁻¹A - c)'x = c_B'B⁻¹b
# =============================================================================

"""
    basic_variable_names(model) -> Vector{String}

Names of the variables whose simplex basis status is `BASIC`.
"""
basic_variable_names(model::Model) =
    [name(v) for v in all_variables(model) if MOI.get(model, MOI.VariableBasisStatus(), v) == MOI.BASIC]

"""
    basic_row_indices(model) -> Vector{Int}

Indices of the (equality) rows whose basis status is `BASIC`. This happens on
degenerate vertices; the slack of such a row is put into the basis so that
the basis matrix stays square. Relies on the row names set by
`create_normal_form_model`.
"""
function basic_row_indices(model::Model)
    rows = all_constraints(model; include_variable_in_set_constraints = false)
    return [parse(Int, name(c)) for c in rows if MOI.get(model, MOI.ConstraintBasisStatus(), c) == MOI.BASIC]
end

"""
    tableau_dataframe(body, rhs, basis_names, var_names) -> DataFrame

Assemble a tableau as a `DataFrame` with columns `BasisVars`, one column per
variable, and `rhs`.
"""
function tableau_dataframe(body, rhs, basis_names, var_names)
    df = DataFrame(body, Symbol.(var_names))
    df.rhs = rhs
    insertcols!(df, 1, :BasisVars => basis_names)
    return df
end

"""
    generate_first_simplex_tableau(norm_model) -> DataFrame

Initial tableau of a normal-form model, with the slack variables as basis
(`z - c'x = c_offset`). It is primal feasible only when `b >= 0`.
"""
function generate_first_simplex_tableau(norm_model::Model)
    data = lp_matrix_data(norm_model)
    A = Matrix(data.A)
    m, ncols = size(A)
    vars = name.(data.variables)
    slacks = vars[ncols-m+1:end]

    body = vcat(A, -data.c')
    rhs = vcat(data.b_upper, data.c_offset)
    return tableau_dataframe(body, rhs, vcat(slacks, "z"), vars)
end

"""
    generate_optimal_simplex_tableau(norm_model)

Build the optimal tableau of a **solved** normal-form model from the basis
reported by the solver. Returns a named tuple with

* `tableau`    – full optimal tableau (`DataFrame`), last row is `z`;
* `basis`      – names of the basic variables (tableau row order);
* `nonbasic`   – names of the non-basic variables;
* `variables`  – names of all variables (structural first, then slacks);
* `structural` – names of the structural (non-slack) variables.
"""
function generate_optimal_simplex_tableau(norm_model::Model)
    data = lp_matrix_data(norm_model)
    A = Matrix(data.A)
    m, ncols = size(A)
    vars = name.(data.variables)
    structural = vars[1:ncols-m]

    basis = complete_basis(A, vars, basic_variable_names(norm_model),
                           vars[length(structural) .+ basic_row_indices(norm_model)])
    basis_idx = indexin(basis, vars)
    c_B = data.c[basis_idx]

    B = lu(A[:, basis_idx])          # factorise once instead of forming B⁻¹
    body = B \ A
    rhs = B \ data.b_upper
    z_row = (c_B' * body)' .- data.c
    z_val = c_B' * rhs + data.c_offset

    tableau = tableau_dataframe(vcat(body, z_row'), vcat(rhs, z_val), vcat(basis, "z"), vars)
    return (tableau = tableau, basis = basis, nonbasic = setdiff(vars, basis),
            variables = vars, structural = structural)
end

"""
    complete_basis(A, vars, basic_vars, basic_row_slacks) -> Vector{String}

Square basis for `[A I]` from the solver's report: its basic variables, plus,
for every row the solver marks basic (degenerate vertex), that row's slack.
If a slack is already basic, or the columns are still too few, the basis is
completed greedily with further slack columns that keep it non-singular.
"""
function complete_basis(A, vars, basic_vars, basic_row_slacks)
    m = size(A, 1)
    basis = unique(vcat(basic_vars, basic_row_slacks))
    rank(A[:, indexin(basis, vars)]) == length(basis) ||
        error("The basis reported by the solver is singular.")
    for v in vars[end-m+1:end]                       # candidate slack columns
        length(basis) == m && break
        v in basis && continue
        candidate = vcat(basis, v)
        rank(A[:, indexin(candidate, vars)]) == length(candidate) && (basis = candidate)
    end
    length(basis) == m || error("Could not recover a square basis ($(length(basis)) columns for $m rows).")
    return basis
end

"""
    verify_optimal_tableau(norm_model; atol = 1e-7)

Check the optimal tableau of a **solved** normal-form model against the
optimality conditions and against the solver itself. Returns a vector of
named tuples `(check, ok, detail)`:

* primal feasibility (`B⁻¹b >= 0`) and dual feasibility (reduced costs `>= 0`);
* the basic values and the objective agree with the solver's solution;
* the simplex multipliers `y = c_B'B⁻¹` agree with the solver's shadow prices;
* basic columns form an identity in the tableau.
"""
function verify_optimal_tableau(norm_model::Model; atol::Real = 1e-7)
    opt = generate_optimal_simplex_tableau(norm_model)
    T = opt.tableau
    m = length(opt.basis)
    data = lp_matrix_data(norm_model)
    A = Matrix(data.A)
    vars = opt.variables
    basis_idx = indexin(opt.basis, vars)
    rhs = T.rhs[1:m]
    d = Vector(T[m+1, vars])                                 # reduced costs
    body = Matrix(T[1:m, vars])

    x_tab = zeros(length(vars)); x_tab[basis_idx] = rhs
    x_sol = value.(all_variables(norm_model))
    y_tab = (lu(A[:, basis_idx])' \ data.c[basis_idx])         # y' = c_B' B⁻¹
    rows = all_constraints(norm_model; include_variable_in_set_constraints = false)
    y_sol = shadow_price.(rows)

    maxerr(a, b) = maximum(abs.(a .- b); init = 0.0)
    fmt(x) = string(round(x; sigdigits = 3))
    return [
        (check = "Primal feasible: every basic value B⁻¹b ≥ 0", ok = minimum(rhs) >= -atol,
         detail = "smallest basic value = $(fmt(minimum(rhs)))"),
        (check = "Optimal: every reduced cost ≥ 0 (maximisation)", ok = minimum(d) >= -atol,
         detail = "smallest reduced cost = $(fmt(minimum(d)))"),
        (check = "Basic columns form an identity matrix", ok = maxerr(body[:, basis_idx], I(m)) < atol,
         detail = "max deviation = $(fmt(maxerr(body[:, basis_idx], I(m))))"),
        (check = "Tableau values match the solver's solution", ok = maxerr(x_tab, x_sol) < 1e-6,
         detail = "max difference = $(fmt(maxerr(x_tab, x_sol)))"),
        (check = "Tableau objective matches the solver's objective", ok = abs(T.rhs[m+1] - objective_value(norm_model)) < 1e-6,
         detail = "tableau $(fmt(T.rhs[m+1])), solver $(fmt(objective_value(norm_model)))"),
        (check = "Simplex multipliers c_B'B⁻¹ match the solver's shadow prices", ok = maxerr(y_tab, y_sol) < 1e-6,
         detail = "max difference = $(fmt(maxerr(y_tab, y_sol)))"),
    ]
end

"""
    compute_rays(tableau, structural, nonbasic) -> Matrix{Float64}

Extreme rays of the LP cone at the optimal vertex, restricted to the
structural variables. Column `j` is the direction obtained by increasing the
`j`-th non-basic variable:

* a basic structural variable moves by `-(B⁻¹N)[row, j]`;
* a non-basic structural variable moves by `1` along its own ray, `0` otherwise.

Rows follow the order of `structural`.
"""
function compute_rays(tableau::DataFrame, structural, nonbasic)
    rays = zeros(length(structural), length(nonbasic))
    basic_rows = Dict(v => i for (i, v) in enumerate(tableau.BasisVars))
    for (i, v) in enumerate(structural)
        if haskey(basic_rows, v)
            rays[i, :] = -Vector(tableau[basic_rows[v], nonbasic])
        else
            rays[i, findfirst(==(v), nonbasic)] = 1.0
        end
    end
    return rays
end
