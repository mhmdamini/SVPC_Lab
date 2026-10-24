# =============================================================================
# One entry point for every cut family, used by the web interface.
# =============================================================================

"""Cut families understood by `generate_cut`, with display names."""
const CUT_METHODS = (
    gfc      = "Gomory fractional cut (GFC)",
    sic      = "Standard intersection cut (SIC)",
    gmi      = "Gomory mixed-integer (GMI)",
    vpc      = "V-polyhedral cut (VPC)",
    vpc_mono = "VPC + monoidal strengthening",
)

"""
    generate_cut(relaxed_model, var_index, integer_vars, method;
                 normalization = :auto, optimizer = HiGHS.Optimizer)
        -> (α, β, note)

Generate one cut `α'x >= β` for the **solved** LP relaxation `relaxed_model`,
using variable `var_index` (the split variable for VPCs and SIC, the tableau row for
GFC and GMI). `method` is one of `keys(CUT_METHODS)`; `normalization` is passed to
`solve_prlp`. `note` is a human-readable remark about how the cut was obtained.
"""
function generate_cut(relaxed_model::Model, var_index::Int, integer_vars, method::Symbol;
                      normalization::Symbol = :auto, optimizer = HiGHS.Optimizer)
    haskey(CUT_METHODS, method) || error("Unknown cut method :$method.")
    std_model = create_standard_form_model(relaxed_model)

    if method == :sic
        cut = sic_cut(std_model, var_index; optimizer)
        return (α = cut.α, β = cut.β, note = join(cut.derivation, "\n"))
    end
    if method == :gfc
        cut = gomory_fractional_cut(std_model, var_index, integer_vars; optimizer)
        return (α = cut.α, β = cut.β, note = join(cut.derivation, "\n"))
    end
    if method == :gmi
        cut = gmi_cut(std_model, var_index, integer_vars; optimizer)
        return (α = cut.α, β = cut.β, note = join(cut.derivation, "\n"))
    end

    # ---- VPC: split disjunction on x_k ----
    subproblems, D0 = create_disjunctive_subproblems(std_model, var_index, integer_vars; optimizer)
    for (t, sp) in enumerate(subproblems)
        optimize_model!(sp; optimizer)
        if termination_status(sp) == MOI.INFEASIBLE
            # Only the other side survives, so its split row is the cut:
            # side 1 empty → x_k >= ⌈x̄⌉;  side 2 empty → -x_k >= -⌊x̄⌋.
            e = zeros(length(all_variables(relaxed_model)))
            e[var_index] = 1.0
            α, β = t == 1 ? (e, -D0[2]) : (-e, -D0[1])
            return (α = α, β = β, note = "Side $t of the split is infeasible, so the other side's split row is the cut.")
        end
    end

    # Points and rays of both sides → point-ray LP.
    results = [process_subproblem(sp; optimizer) for sp in subproblems]
    ω, _ = get_variable_values(relaxed_model)
    prlp, α, β = solve_prlp(ω, [r.point for r in results], [r.rays for r in results];
                            normalization, optimizer, verbose = false)
    note = is_fixed(prlp[:β]) ? "PRLP solved with β = 1." :
                                "PRLP solved with ‖α‖₁ ≤ 1 (the β = 1 version had no optimal solution)."

    if method == :vpc_mono
        U = prepare_monoidal_data(subproblems, results, α)
        α̃ = monoidal_strengthening(subproblems, var_index, U, α, integer_vars)
        note *= α̃ ≈ α ? " Monoidal strengthening left the coefficients unchanged." :
                         " Coefficients strengthened monoidally."
        α = α̃
    end
    return (α = α, β = β, note = note)
end
