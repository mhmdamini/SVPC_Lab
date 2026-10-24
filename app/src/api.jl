# =============================================================================
# The JSON API behind the interface: parse a problem, generate the requested
# cuts, and package bounds, solutions, the tableau and plot geometry.
# =============================================================================

# ---- Parsing ----------------------------------------------------------------

"""Build a JuMP model from LP- or MPS-format text."""
function parse_model(text::AbstractString, format::AbstractString)
    fmt = lowercase(format) == "mps" ? MOI.FileFormats.FORMAT_MPS : MOI.FileFormats.FORMAT_LP
    src = MOI.FileFormats.Model(format = fmt)
    read!(IOBuffer(text), src)
    model = Model()
    MOI.copy_to(backend(model), src)
    return model
end

# ---- Formatting helpers -----------------------------------------------------

"""Round for display and turn non-finite numbers into `nothing` (valid JSON)."""
clean(x::Real) = isfinite(x) ? round(Float64(x); sigdigits = 6) + 0.0 : nothing
clean(x::AbstractVector) = [clean(v) for v in x]
clean(x::AbstractMatrix) = [clean(x[i, :]) for i in axes(x, 1)]
clean(x) = x

"""Pretty-print `α'x >= β` with variable names, e.g. `-0.2 x1 + 0.8 x2 ≥ 0.2`."""
function format_cut(α, β, names)
    terms = String[]
    for (a, nm) in zip(α, names)
        abs(a) < 1e-9 && continue
        c = round(a; sigdigits = 4)
        mag = abs(c) == 1 ? "" : string(isinteger(c) ? Int(abs(c)) : abs(c), " ")
        push!(terms, (isempty(terms) ? (c < 0 ? "−" : "") : (c < 0 ? " − " : " + ")) * mag * nm)
    end
    lhs = isempty(terms) ? "0" : join(terms)
    rhs = round(β; sigdigits = 4)
    rhs_text = replace(string(isinteger(rhs) ? Int(rhs) : rhs), "-" => "−")
    return "$lhs ≥ $rhs_text"
end

# ---- Simplex report -----------------------------------------------------------

# Iterations are only sent for problems small enough to read as tables.
const MAX_TABLEAU_ROWS = 40
const MAX_TABLEAU_COLS = 60

"""JSON-ready view of one recorded simplex step."""
function step_json(st::SimplexStep)
    return Dict(
        "phase" => st.phase, "columns" => st.columns, "basis" => st.basis,
        "rows" => clean(st.body), "rhs" => clean(st.rhs), "reduced" => clean(st.reduced),
        "objective" => clean(st.objective),
        "entering" => st.entering === nothing ? nothing : st.entering - 1,   # 0-based for JS
        "leaving" => st.leaving === nothing ? nothing : st.leaving - 1,
        "ratios" => [r === nothing ? nothing : clean(r) for r in st.ratios],
        "note" => st.note)
end

"""
    simplex_report(relaxed_model, sense) -> Dict

Simplex view of the LP relaxation in standard form (max c'x, Ax <= b, x >= 0):
textbook two-phase iterations, and an analysis of the optimal basis computed
from the solver's basis (checks, values, reduced costs, shadow prices,
degeneracy, alternative optima). Slacks are named `s1, s2, …` after their rows.
"""
function simplex_report(relaxed_model, sense)
    std = create_standard_form_model(relaxed_model)
    sdata = lp_matrix_data(std)
    m, n = size(sdata.A)
    row_labels = [name(c) for c in sdata.affine_constraints]
    var_names = name.(sdata.variables)

    # Textbook two-phase simplex (independent of the solver).
    run = tableau_simplex(std)
    too_big = m > MAX_TABLEAU_ROWS || n + m > MAX_TABLEAU_COLS
    iterations = too_big ? nothing : [step_json(st) for st in run.steps]

    # Optimal tableau from the solver's basis, and its verification.
    normal = create_normal_form_model(std)
    optimize_model!(normal)
    opt = generate_optimal_simplex_tableau(normal)
    pretty(v) = startswith(v, "β[") ? var_names[parse(Int, v[4:end-1])] :
                 startswith(v, "s[") ? "s" * v[3:end-1] : v
    T = opt.tableau
    vars = opt.variables
    nb = length(opt.basis)
    values = zeros(length(vars))
    for (r, v) in enumerate(opt.basis)
        values[findfirst(==(v), vars)] = T.rhs[r]
    end
    reduced = Vector(T[nb+1, vars])
    rows = all_constraints(normal; include_variable_in_set_constraints = false)
    duals = shadow_price.(rows)

    checks = [Dict("check" => c.check, "ok" => c.ok, "detail" => c.detail) for c in verify_optimal_tableau(normal)]
    push!(checks, Dict(
        "check" => "Textbook two-phase simplex reaches the same optimum",
        "ok" => run.status == :optimal && isapprox(run.objective, objective_value(normal); atol = 1e-6),
        "detail" => run.status == :optimal ?
            "$(count(s -> s.entering !== nothing, run.steps)) pivots ($(count(s -> s.phase == 1 && s.entering !== nothing, run.steps)) in Phase I); objective $(round(run.objective; sigdigits = 6))" :
            "status: $(run.status)"))

    tol = 1e-9
    variables = [Dict("name" => pretty(v), "kind" => startswith(v, "s[") ? "slack" : "structural",
                      "basic" => v in opt.basis, "value" => clean(values[j]), "reduced_cost" => clean(reduced[j]))
                 for (j, v) in enumerate(vars)]
    constraints = [Dict("row" => i, "label" => row_labels[i], "slack" => "s$i",
                        "slack_value" => clean(values[n+i]), "binding" => abs(values[n+i]) < 1e-7,
                        "shadow_price" => clean(duals[i])) for i in 1:m]

    return Dict(
        "status" => String(run.status),
        "sense_note" => sense == "min" ?
            "The original problem minimises, so the standard form maximises −(objective); z in the tableau is −(original objective)." :
            "The standard form maximises the original objective.",
        "row_labels" => row_labels,
        "iterations" => iterations,
        "too_big" => too_big,
        "optimal" => Dict("columns" => pretty.(vars), "basis" => pretty.(opt.basis),
                          "rows" => clean(Matrix(T[1:nb, vars])), "rhs" => clean(T.rhs[1:nb]),
                          "reduced" => clean(reduced), "objective" => clean(T.rhs[nb+1])),
        "checks" => checks,
        "variables" => variables,
        "constraints" => constraints,
        "degenerate" => [pretty(v) for (r, v) in enumerate(opt.basis) if abs(T.rhs[r]) < 1e-7],
        "alternative" => [pretty(v) for (j, v) in enumerate(vars) if !(v in opt.basis) && abs(reduced[j]) < 1e-7],
    )
end

# ---- Main handler -----------------------------------------------------------

function solve_request(req)
    model = parse_model(req.text, get(req, :format, "lp"))
    vars = all_variables(model)
    names_ = name.(vars)
    n = length(vars)
    data = lp_matrix_data(model)
    integer_vars = integer_variable_indices(model)
    isempty(integer_vars) && error("The problem has no integer variables, so there is nothing to cut.")
    any(<(0), data.x_lower) &&
        error("All variables must have lower bound ≥ 0 (the cut generator works in standard form x ≥ 0).")

    # MIP optimum (for the gap) and LP relaxation.
    optimize_model!(model)
    mip_ok = termination_status(model) == MOI.OPTIMAL
    z_int = mip_ok ? objective_value(model) : nothing
    x_int = mip_ok ? value.(vars) : nothing

    relaxed = relax_model(model)
    optimize_model!(relaxed)
    termination_status(relaxed) == MOI.OPTIMAL ||
        error("The LP relaxation has no optimal solution (status: $(termination_status(relaxed))).")
    z_lp = objective_value(relaxed)
    ω, _ = get_variable_values(relaxed)

    fractional = fractional_integer_variables(relaxed, integer_vars)
    requested = get(req, :var_index, nothing)
    k = requested === nothing ? (isempty(fractional) ? nothing : first(fractional)) : Int(requested)
    k === nothing || k in fractional ||
        error("Variable $(names_[k]) is not fractional in the LP optimum; choose one of the listed variables.")

    # Cuts.
    normalization = Symbol(get(req, :normalization, "auto"))
    cuts = Dict{String,Any}[]
    good_cuts = []
    if k !== nothing
        for m in get(req, :cuts, ["vpc"])
            method = Symbol(m)
            entry = Dict{String,Any}("method" => String(method), "label" => CUT_METHODS[method])
            try
                cut = generate_cut(relaxed, k, integer_vars, method; normalization)
                after = add_cut(relaxed, cut.α, cut.β; verbose = false)
                optimize_model!(after)
                z_cut = objective_value(after)
                merge!(entry, Dict(
                    "alpha" => clean(cut.α), "beta" => clean(cut.β),
                    "text" => format_cut(cut.α, cut.β, names_), "note" => cut.note,
                    "separates" => cut.α' * ω < cut.β - 1e-6,
                    "bound" => clean(z_cut), "solution" => clean(value.(all_variables(after))),
                    "gap_closed" => z_int === nothing ? nothing : clean(gap_closed(z_cut, z_int, z_lp))))
                push!(good_cuts, (α = cut.α, β = cut.β, entry = entry))
            catch err
                entry["error"] = sprint(showerror, err)
            end
            push!(cuts, entry)
        end
    end

    # Geometry (2-D / 3-D only) and validity check on the enumerated integer points.
    geometry = nothing
    if n in (2, 3)
        all_int = length(integer_vars) == n
        geometry = region_geometry(relaxed, ω, good_cuts; all_integer = all_int)
        for c in good_cuts
            pts = geometry["integer_feasible"]
            isempty(pts) && continue
            c.entry["violated_points"] = count(p -> c.α' * p < c.β - 1e-6, pts)
            c.entry["checked_points"] = length(pts)
        end
        geometry = Dict(key => clean_any(v) for (key, v) in geometry)
    end

    return Dict(
        "variables" => names_, "sense" => data.sense == MOI.MIN_SENSE ? "min" : "max",
        "model_text" => sprint(print, model),
        "integer_vars" => integer_vars,
        "lp" => Dict("bound" => clean(z_lp), "solution" => clean(ω)),
        "mip" => Dict("status" => string(termination_status(model)), "bound" => clean_any(z_int),
                      "solution" => clean_any(x_int)),
        "fractional" => [Dict("index" => j, "name" => names_[j], "value" => clean(ω[j])) for j in fractional],
        "split_var" => k,
        "cuts" => cuts,
        "simplex" => simplex_report(relaxed, data.sense == MOI.MIN_SENSE ? "min" : "max"),
        "geometry" => geometry,
    )
end

clean_any(x::Nothing) = nothing
clean_any(x::Real) = clean(x)
clean_any(x::AbstractVector{<:Real}) = clean(x)
clean_any(x::AbstractVector) = [clean_any(v) for v in x]
clean_any(x) = x
