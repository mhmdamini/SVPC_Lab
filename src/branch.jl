# =============================================================================
# Branch and bound, and branch and cut.
#
# Each node is the LP relaxation with extra variable bounds from branching.
# A node is pruned when its LP is infeasible, when its bound cannot beat the
# incumbent, or when its LP solution is integral (new incumbent). Otherwise the
# branching variable x_k (fractional value x̄_k) creates two children
#
#     x_k <= ⌊x̄_k⌋        and        x_k >= ⌈x̄_k⌉.
#
# Branch and cut adds cutting-plane rounds (`separation_round!`):
#   * at the root — those cuts are valid globally and are added to every node;
#   * optionally at every node — those cuts rely on the node's branching
#     bounds, so they are kept only for that node's subtree (local cuts).
# =============================================================================

"""One node of the branch-and-bound tree."""
mutable struct BBNode
    id::Int
    parent::Int                        # 0 for the root
    depth::Int
    branch::String                     # e.g. "x1 ≤ 0" (empty for the root)
    lower::Vector{Float64}             # variable bounds at this node
    upper::Vector{Float64}
    parent_bound::Float64              # LP bound inherited from the parent
    bound::Union{Nothing,Float64}      # LP bound of this node (after its cuts)
    x::Union{Nothing,Vector{Float64}}
    status::Symbol                     # :open, :branched, :integer, :infeasible, :pruned_bound
    local_cuts::Vector{AddedCut}       # cuts valid in this subtree only
    cuts_added::Int                    # cuts generated at this node
    order::Int                         # processing order (0 = not processed)
end

"""
    branch_and_bound(model; node_selection = :best_bound, branching = :most_fractional,
                     cut_methods = Symbol[], root_cut_rounds = 0, node_cut_rounds = 0,
                     split_rule = :first, max_nodes = 1000, optimizer = HiGHS.Optimizer)

Solve the MIP `model` (not modified) by LP-based branch and bound. With
`cut_methods` and `root_cut_rounds > 0` and/or `node_cut_rounds > 0` it is
branch and cut. Returns a named tuple

* `status` – `:optimal`, `:infeasible`, `:node_limit` or `:root_lp_failed`;
* `objective`, `x` – best integer solution found (`nothing` if none);
* `nodes` – every `BBNode` (ids 1, 2, …; root is 1);
* `history` – after each processed node: `(node, incumbent, best_bound)`;
* `root_bound` (before cuts) and `root_bound_after_cuts`, `global_cuts`.

Options: `node_selection` ∈ `:best_bound, :depth_first, :breadth_first`;
`branching` ∈ `:most_fractional, :first_fractional`.
"""
function branch_and_bound(model::Model; node_selection::Symbol = :best_bound, branching::Symbol = :most_fractional,
                          cut_methods = Symbol[], root_cut_rounds::Int = 0, node_cut_rounds::Int = 0,
                          split_rule::Symbol = :first, normalization::Symbol = :auto, max_nodes::Int = 1000,
                          optimizer = HiGHS.Optimizer, tol::Real = 1e-6)
    node_selection in (:best_bound, :depth_first, :breadth_first) ||
        error("Unknown node selection :$node_selection.")
    branching in (:most_fractional, :first_fractional) || error("Unknown branching rule :$branching.")

    base = relax_model(model)
    data = lp_matrix_data(base)
    vars = all_variables(base)
    names_ = name.(vars)
    n = length(vars)
    integer_vars = integer_variable_indices(model)
    minimize = data.sense == MOI.MIN_SENSE
    better(a, b) = minimize ? a < b - tol : a > b + tol        # a strictly better than b
    worst = minimize ? Inf : -Inf

    global_cuts = AddedCut[]
    nodes = BBNode[]
    history = NamedTuple[]
    incumbent, incumbent_x = nothing, nothing

    # LP of a node: base relaxation + node bounds + global cuts + the node's local cuts.
    function node_lp(node)
        lp = copy(base)
        v = all_variables(lp)
        for j in 1:n
            isfinite(node.lower[j]) ? set_lower_bound(v[j], node.lower[j]) : (has_lower_bound(v[j]) && delete_lower_bound(v[j]))
            isfinite(node.upper[j]) ? set_upper_bound(v[j], node.upper[j]) : (has_upper_bound(v[j]) && delete_upper_bound(v[j]))
        end
        for c in Iterators.flatten((global_cuts, node.local_cuts))
            @constraint(lp, c.α' * v >= c.β)
        end
        return lp
    end

    # Cut rounds on a solved node LP; returns the new cuts (lp re-solved after each round).
    function cut_rounds!(lp, rounds, label)
        found = AddedCut[]
        for _ in 1:rounds
            isempty(fractional_integer_variables(lp, integer_vars; atol = tol)) && break
            cuts, _ = separation_round!(lp, integer_vars, cut_methods; split_rule, label, normalization, optimizer, tol)
            isempty(cuts) && break
            append!(found, cuts)
            optimize_model!(lp; optimizer)
            termination_status(lp) == MOI.OPTIMAL || break
        end
        return found
    end

    branch_var(x) = begin
        frac = fractional_integer_variables_at(x, integer_vars, tol)
        isempty(frac) ? nothing :
        branching == :first_fractional ? first(frac) :
        frac[argmin([abs(x[j] - floor(x[j]) - 0.5) for j in frac])]
    end

    best_open_bound(open) = isempty(open) ? (incumbent === nothing ? worst : incumbent) :
        (minimize ? minimum(nd.parent_bound for nd in open) : maximum(nd.parent_bound for nd in open))

    # ---- Root ----
    root = BBNode(1, 0, 0, "", [isfinite(l) ? l : -Inf for l in data.x_lower], [isfinite(u) ? u : Inf for u in data.x_upper],
                  worst == Inf ? -Inf : Inf, nothing, nothing, :open, AddedCut[], 0, 0)
    push!(nodes, root)
    open = BBNode[root]
    root_bound, root_bound_after_cuts = nothing, nothing
    processed = 0

    while !isempty(open)
        processed >= max_nodes && break
        # Node selection.
        idx = node_selection == :depth_first ? length(open) :
              node_selection == :breadth_first ? 1 :
              (minimize ? argmin([nd.parent_bound for nd in open]) : argmax([nd.parent_bound for nd in open]))
        node = open[idx]
        deleteat!(open, idx)

        # Prune by bound before solving (the incumbent may have improved since it was created).
        if incumbent !== nothing && node.id != 1 && !better(node.parent_bound, incumbent)
            node.status = :pruned_bound
            continue
        end

        processed += 1
        node.order = processed
        lp = node_lp(node)
        optimize_model!(lp; optimizer)
        st = termination_status(lp)
        if st in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
            node.status = :infeasible
        elseif st != MOI.OPTIMAL
            node.id == 1 && return (status = :root_lp_failed, objective = nothing, x = nothing, nodes = nodes,
                                    history = history, root_bound = nothing, root_bound_after_cuts = nothing,
                                    global_cuts = global_cuts)
            node.status = :infeasible
        else
            node.id == 1 && (root_bound = objective_value(lp))
            # Cutting planes: global at the root, local below it.
            rounds = node.id == 1 ? root_cut_rounds : node_cut_rounds
            if !isempty(cut_methods) && rounds > 0
                new_cuts = cut_rounds!(lp, rounds, node.id)
                node.cuts_added = length(new_cuts)
                node.id == 1 ? append!(global_cuts, new_cuts) : append!(node.local_cuts, new_cuts)
                if termination_status(lp) != MOI.OPTIMAL
                    optimize_model!(lp; optimizer)
                end
            end
            if termination_status(lp) in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
                node.status = :infeasible
            else
                node.bound = objective_value(lp)
                node.x = value.(all_variables(lp))
                node.id == 1 && (root_bound_after_cuts = node.bound)
                k = branch_var(node.x)
                if incumbent !== nothing && !better(node.bound, incumbent)
                    node.status = :pruned_bound
                elseif k === nothing
                    node.status = :integer
                    if incumbent === nothing || better(node.bound, incumbent)
                        incumbent, incumbent_x = node.bound, round.(node.x; digits = 9)
                    end
                else
                    node.status = :branched
                    xk = node.x[k]
                    for (lo, up, text) in ((node.lower[k], floor(xk), "$(names_[k]) ≤ $(number_text(floor(xk)))"),
                                           (ceil(xk), node.upper[k], "$(names_[k]) ≥ $(number_text(ceil(xk)))"))
                        child_lower, child_upper = copy(node.lower), copy(node.upper)
                        child_lower[k], child_upper[k] = lo, up
                        child = BBNode(length(nodes) + 1, node.id, node.depth + 1, text, child_lower, child_upper,
                                       node.bound, nothing, nothing, :open, copy(node.local_cuts), 0, 0)
                        push!(nodes, child)
                        lo > up ? (child.status = :infeasible) : push!(open, child)
                    end
                end
            end
        end
        push!(history, (node = processed, incumbent = incumbent, best_bound = best_open_bound(open)))
    end

    status = isempty(open) ? (incumbent === nothing ? :infeasible : :optimal) : :node_limit
    return (status = status, objective = incumbent, x = incumbent_x, nodes = nodes, history = history,
            root_bound = root_bound, root_bound_after_cuts = root_bound_after_cuts, global_cuts = global_cuts)
end

"""Indices in `integer_vars` whose value in `x` is fractional (beyond `tol`)."""
fractional_integer_variables_at(x, integer_vars, tol) =
    [j for j in integer_vars if abs(x[j] - round(x[j])) > tol]
