# =============================================================================
# Iterative cutting-plane algorithm.
#
# Repeat until the LP solution is integral (then it is optimal for the MIP),
# no cut separates it any more, the bound stalls, or a round limit is hit:
#
#   1. solve the current relaxation (LP + all cuts added so far);
#   2. pick split variables among the fractional integer variables;
#   3. generate one cut per (variable, cut family) with `generate_cut`;
#   4. add every cut that separates the LP solution.
#
# With one family this is e.g. Gomory's classical cutting-plane method (GFC
# rounds) or a GMI loop; with several families their cuts are added together.
#
# For pure-integer problems each new cut is scaled to integer coefficients and
# its right-hand side rounded up (Chvátal–Gomory rounding): valid because α'x is
# an integer at integer points, slightly stronger, and it keeps the data
# integral so Gomory fractional cuts remain exact in later rounds.
# =============================================================================

"""One cut added by an iterative method."""
struct AddedCut
    round::Int              # cutting-plane round (or B&B node id) that produced it
    method::Symbol          # cut family (`keys(CUT_METHODS)`)
    var_index::Int          # split variable / tableau row it came from
    α::Vector{Float64}
    β::Float64
end

"""
    integer_cut(α, β; max_denominator = 1_000_000) -> (α, β, ok)

Scale `α'x >= β` to integer coefficients and round `β` up (Chvátal–Gomory).
Only valid when every variable is integer. Returns `ok = false` (and the cut
unchanged) when the coefficients are not rational with a small denominator.
"""
function integer_cut(α::AbstractVector, β::Real; max_denominator::Int = 1_000_000)
    scale = 1
    for a in α
        abs(a) < 1e-12 && continue
        r = rationalize(Int, a / maximum(abs, α); tol = 1e-10)
        denominator(r) > max_denominator && return (collect(float.(α)), float(β), false)
        scale = lcm(scale, denominator(r))
        scale > max_denominator && return (collect(float.(α)), float(β), false)
    end
    s = scale / maximum(abs, α)
    α_int = round.(α .* s)
    maximum(abs.(α_int .- α .* s)) < 1e-6 || return (collect(float.(α)), float(β), false)
    return (α_int, ceil(β * s - 1e-9), true)
end

"""
    choose_split_variables(x, fractional, rule) -> Vector{Int}

Which fractional variables to cut on in a round: `:first`, `:most_fractional`
(value closest to .5) or `:all`.
"""
function choose_split_variables(x, fractional, rule::Symbol)
    isempty(fractional) && return Int[]
    rule == :all && return collect(fractional)
    rule == :first && return [first(fractional)]
    rule == :most_fractional && return [fractional[argmin([abs(x[j] - floor(x[j]) - 0.5) for j in fractional])]]
    error("Unknown split rule :$rule (use :first, :most_fractional or :all).")
end

"""
    separation_round!(lp, integer_vars, methods; split_rule, label, normalization, optimizer, tol)
        -> (cuts::Vector{AddedCut}, errors::Vector{String})

One round of cut generation on the **solved** LP `lp`: generate a cut for
every chosen split variable and every family in `methods`, keep those that
separate the LP solution, and add them to `lp` (not re-solved). `label` is
stored as the cut's round number.
"""
function separation_round!(lp::Model, integer_vars, methods; split_rule::Symbol = :first, label::Int = 0,
                           normalization::Symbol = :auto, optimizer = HiGHS.Optimizer, tol::Real = 1e-6)
    x, names_ = get_variable_values(lp)
    vars = all_variables(lp)
    pure_integer = length(integer_vars) == length(vars)
    fractional = fractional_integer_variables(lp, integer_vars; atol = tol)
    cuts, errors = AddedCut[], String[]
    for k in choose_split_variables(x, fractional, split_rule), m in methods
        try
            cut = generate_cut(lp, k, integer_vars, m; normalization, optimizer)
            α, β = cut.α, cut.β
            if pure_integer
                αi, βi, ok = integer_cut(α, β)
                ok && ((α, β) = (αi, βi))
            end
            if α' * x < β - tol * max(1.0, maximum(abs, α))       # keep only separating cuts
                push!(cuts, AddedCut(label, m, k, α, β))
            else
                push!(errors, "$(CUT_METHODS[m]) on $(names_[k]): does not cut off the LP solution")
            end
        catch err
            push!(errors, "$(CUT_METHODS[m]) on $(names_[k]): " * first(split(sprint(showerror, err), '\n')))
        end
    end
    for (i, c) in enumerate(cuts)
        con = @constraint(lp, c.α' * vars >= c.β)
        set_name(con, "cut$(label)_$(c.method)_$(i)")
    end
    return cuts, errors
end

"""
    cutting_plane_loop(model; methods = [:gmi], max_rounds = 20, split_rule = :first,
                       stall_rounds = 5, normalization = :auto, optimizer = HiGHS.Optimizer)

Iterative cutting-plane method on the MIP `model` (not modified). Each round
solves the relaxation with every cut so far, then adds the separating cuts of
the families in `methods` (one family: a pure cutting-plane method; several:
their cuts together). Returns a named tuple

* `status` – `:optimal` (LP solution integral ⇒ MIP optimum), `:no_cuts`
  (nothing separates any more), `:stalled` (no bound progress for
  `stall_rounds` rounds), `:round_limit` or `:lp_failed`;
* `rounds` – one entry per round: `round`, `bound`, `x`, `cuts` (added that
  round), `errors` (skipped cuts and why);
* `cuts` – every `AddedCut`; `lp` – the final relaxation with all cuts;
* `objective`, `x`, `integral` – the final LP solution.
"""
function cutting_plane_loop(model::Model; methods = [:gmi], max_rounds::Int = 20, split_rule::Symbol = :first,
                            stall_rounds::Int = 5, normalization::Symbol = :auto,
                            optimizer = HiGHS.Optimizer, tol::Real = 1e-6)
    isempty(methods) && error("Select at least one cut family.")
    integer_vars = integer_variable_indices(model)
    lp = relax_model(model)
    optimize_model!(lp; optimizer)
    termination_status(lp) == MOI.OPTIMAL || error("The LP relaxation has no optimal solution ($(termination_status(lp))).")

    rounds = Any[]
    all_cuts = AddedCut[]
    record(r, cuts, errors) = push!(rounds, (round = r, bound = objective_value(lp), x = value.(all_variables(lp)),
                                             cuts = cuts, errors = errors))
    record(0, AddedCut[], String[])

    status, stall = :round_limit, 0
    for r in 1:max_rounds
        if isempty(fractional_integer_variables(lp, integer_vars; atol = tol))
            status = :optimal
            break
        end
        cuts, errors = separation_round!(lp, integer_vars, methods; split_rule, label = r, normalization, optimizer, tol)
        if isempty(cuts)
            rounds[end] = merge(rounds[end], (errors = vcat(rounds[end].errors, errors),))
            status = :no_cuts
            break
        end
        previous = rounds[end].bound                     # (lp was modified, so read the record)
        optimize_model!(lp; optimizer)
        if termination_status(lp) != MOI.OPTIMAL
            status = :lp_failed
            break
        end
        append!(all_cuts, cuts)
        record(r, cuts, errors)
        stall = abs(objective_value(lp) - previous) <= 1e-9 * max(1, abs(previous)) ? stall + 1 : 0
        if stall >= stall_rounds
            status = :stalled
            break
        end
    end
    # The last recorded round is the last relaxation that was solved to optimality.
    last_round = rounds[end]
    integral = isempty(intersect(find_fractional_elements(last_round.x; atol = tol)[2], integer_vars))
    status == :round_limit && integral && (status = :optimal)
    return (status = status, rounds = rounds, cuts = all_cuts, lp = lp, objective = last_round.bound,
            x = last_round.x, integral = integral)
end
