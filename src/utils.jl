# =============================================================================
# Basic model utilities: reading, relaxing, solving and reporting.
# =============================================================================

"""
    read_model(path) -> Model

Read a model from any file format JuMP understands (`.mps`, `.mps.gz`, `.lp`, ...).
"""
read_model(path::AbstractString) = read_from_file(path)

"""
    integer_variable_indices(model) -> Vector{Int}

Indices (in `all_variables` order) of the integer **and** binary variables.
JuMP's `lp_matrix_data` lists the two kinds separately.
"""
function integer_variable_indices(model::Model)
    data = lp_matrix_data(model)
    return sort(union(data.integers, data.binaries))
end

"""
    relax_model(model) -> Model

Return a copy of `model` with every integrality/binary restriction removed,
i.e. its LP relaxation. The original model is left untouched.
"""
function relax_model(model::Model)
    relaxed = copy(model)
    relax_integrality(relaxed)
    return relaxed
end

"""
    optimize_model!(model; optimizer = HiGHS.Optimizer, verbose = false) -> model

Attach `optimizer` to `model` (needed after `copy`, which drops the solver)
and solve it. Solver output is suppressed unless `verbose = true`.
"""
function optimize_model!(model::Model; optimizer = HiGHS.Optimizer, verbose::Bool = false)
    set_optimizer(model, optimizer)
    verbose ? unset_silent(model) : set_silent(model)
    optimize!(model)
    return model
end

"""
    get_variable_values(model) -> (values::Vector{Float64}, names::Vector{String})

Primal values and names of all variables, in `all_variables(model)` order.
Throws if the model was not solved to optimality.
"""
function get_variable_values(model::Model)
    status = termination_status(model)
    status == MOI.OPTIMAL ||
        error("Model was not solved to optimality (termination status: $status).")
    vars = all_variables(model)
    return value.(vars), name.(vars)
end

"""
    retrieve_objective_values(model) -> (objective_value, termination_status)
"""
retrieve_objective_values(model::Model) = (objective_value(model), termination_status(model))

"""
    gap_closed(z_cut, z_int, z_lp) -> Float64

Percentage of the integrality gap closed by a cut:

    100 * (z_cut - z_lp) / (z_int - z_lp)

where `z_lp` is the LP-relaxation bound, `z_int` the optimal MIP value and
`z_cut` the LP bound after adding the cut. Returns 100 when the relaxation is
already tight (`z_int == z_lp`).
"""
function gap_closed(z_cut::Real, z_int::Real, z_lp::Real)
    isapprox(z_int, z_lp; atol = 1e-9) && return 100.0
    return 100 * (z_cut - z_lp) / (z_int - z_lp)
end

"""
    print_variable_values(values, names)

Print one `name = value` line per variable.
"""
function print_variable_values(values, names)
    for (n, v) in zip(names, values)
        println("  $n = $v")
    end
end

"""
    print_banner(message)

Print `message` framed in a box, used to separate the stages of the pipeline.
"""
function print_banner(message::AbstractString)
    width = length(message) + 2
    println()
    println("╔", "═"^width, "╗")
    println("║ ", message, " ║")
    println("╚", "═"^width, "╝")
end
