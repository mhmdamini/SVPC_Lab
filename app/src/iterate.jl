# =============================================================================
# JSON API for the Solve tab (POST /api/iterate): iterative cutting planes,
# branch and bound, and branch and cut on a problem given as LP/MPS text.
# =============================================================================

const MAX_ROUNDS_LIMIT = 200
const MAX_NODES_LIMIT = 5000

"""Common header of an iterate response: variables, LP bound, MIP optimum."""
function problem_summary(model, relaxed)
    vars = all_variables(model)
    data = lp_matrix_data(model)
    mip = copy(model)
    optimize_model!(mip)
    mip_ok = termination_status(mip) == MOI.OPTIMAL
    return Dict(
        "variables" => name.(vars), "sense" => data.sense == MOI.MIN_SENSE ? "min" : "max",
        "integer_vars" => integer_variable_indices(model),
        "lp" => Dict("bound" => clean(objective_value(relaxed)), "solution" => clean(value.(all_variables(relaxed)))),
        "mip" => Dict("status" => string(termination_status(mip)), "bound" => mip_ok ? clean(objective_value(mip)) : nothing,
                      "solution" => mip_ok ? clean(value.(all_variables(mip))) : nothing))
end

"""Cut as JSON: inequality text, family, round and source variable."""
cut_json(c::AddedCut, names_) = Dict("text" => format_cut(c.α, c.β, names_), "method" => String(c.method),
                                     "label" => CUT_METHODS[c.method], "round" => c.round,
                                     "var" => names_[c.var_index])

function iterate_request(req)
    model = parse_model(req.text, get(req, :format, "lp"))
    names_ = name.(all_variables(model))
    integer_vars = integer_variable_indices(model)
    isempty(integer_vars) && error("The problem has no integer variables.")
    any(<(0), lp_matrix_data(model).x_lower) && error("All variables must have lower bound ≥ 0.")

    relaxed = relax_model(model)
    optimize_model!(relaxed)
    termination_status(relaxed) == MOI.OPTIMAL ||
        error("The LP relaxation has no optimal solution (status: $(termination_status(relaxed))).")
    out = problem_summary(model, relaxed)
    z_lp, z_ip = objective_value(relaxed), out["mip"]["bound"]
    gap(z) = z_ip === nothing ? nothing : clean(gap_closed(z, z_ip, z_lp))

    algorithm = String(get(req, :algorithm, "cutting_planes"))
    methods = Symbol.(collect(get(req, :methods, ["gmi"])))
    normalization = Symbol(get(req, :normalization, "auto"))
    split_rule = Symbol(get(req, :split_rule, "first"))
    all_integer = length(integer_vars) == length(names_)
    out["algorithm"] = algorithm
    t0 = time()

    if algorithm == "cutting_planes"
        max_rounds = clamp(Int(get(req, :max_rounds, 20)), 1, MAX_ROUNDS_LIMIT)
        r = cutting_plane_loop(model; methods, max_rounds, split_rule, normalization)
        out["status"] = String(r.status)
        out["objective"] = clean(r.objective)
        out["x"] = clean(r.x)
        out["integral"] = r.integral
        out["rounds"] = [Dict("round" => rd.round, "bound" => clean(rd.bound), "gap_closed" => gap(rd.bound),
                              "x" => clean(rd.x), "cuts" => [cut_json(c, names_) for c in rd.cuts],
                              "errors" => rd.errors) for rd in r.rounds]
        out["total_cuts"] = length(r.cuts)
        if length(names_) in (2, 3)
            geo = iteration_geometry(relaxed, [(α = c.α, β = c.β, round = c.round, method = c.method) for c in r.cuts],
                                     length(r.rounds) - 1; all_integer)
            out["geometry"] = geo === nothing ? nothing : Dict(k => clean_any(v) for (k, v) in geo)
        end
    elseif algorithm in ("branch_and_bound", "branch_and_cut")
        cut = algorithm == "branch_and_cut"
        r = branch_and_bound(model;
                             node_selection = Symbol(get(req, :node_selection, "best_bound")),
                             branching = Symbol(get(req, :branching, "most_fractional")),
                             max_nodes = clamp(Int(get(req, :max_nodes, 500)), 1, MAX_NODES_LIMIT),
                             cut_methods = cut ? methods : Symbol[],
                             root_cut_rounds = cut ? clamp(Int(get(req, :root_cut_rounds, 3)), 0, MAX_ROUNDS_LIMIT) : 0,
                             node_cut_rounds = cut ? clamp(Int(get(req, :node_cut_rounds, 0)), 0, 20) : 0,
                             split_rule, normalization)
        out["status"] = String(r.status)
        out["objective"] = clean_any(r.objective)
        out["x"] = clean_any(r.x)
        out["root_bound"] = clean_any(r.root_bound)
        out["root_bound_after_cuts"] = clean_any(r.root_bound_after_cuts)
        out["root_gap_closed"] = r.root_bound_after_cuts === nothing ? nothing : gap(r.root_bound_after_cuts)
        out["nodes"] = [Dict("id" => nd.id, "parent" => nd.parent, "depth" => nd.depth, "branch" => nd.branch,
                             "bound" => clean_any(nd.bound), "status" => String(nd.status), "order" => nd.order,
                             "x" => clean_any(nd.x), "cuts" => nd.cuts_added) for nd in r.nodes]
        out["history"] = [Dict("node" => h.node, "incumbent" => clean_any(h.incumbent), "best_bound" => clean_any(h.best_bound))
                          for h in r.history]
        out["global_cuts"] = [cut_json(c, names_) for c in r.global_cuts]
        out["total_cuts"] = length(r.global_cuts) + sum(nd.cuts_added for nd in r.nodes if nd.id != 1; init = 0)
        if length(names_) in (2, 3)
            geo = iteration_geometry(relaxed, [(α = c.α, β = c.β, round = 1, method = c.method) for c in r.global_cuts],
                                     1; all_integer)
            out["geometry"] = geo === nothing ? nothing : Dict(k => clean_any(v) for (k, v) in geo)
        end
    else
        error("Unknown algorithm \"$algorithm\".")
    end
    out["seconds"] = round(time() - t0; digits = 3)
    return out
end
