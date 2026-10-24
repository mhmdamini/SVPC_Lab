# =============================================================================
# End-to-end run: generate one V-polyhedral cut (VPC) for a MIP, report how
# much of the integrality gap it closes, then try to strengthen it with
# monoidal strengthening.
#
# Usage (from the repository root):
#   julia --project=. run.jl                       # built-in toy example
#   julia --project=. run.jl model.mps             # split on the first fractional variable
#   julia --project=. run.jl model.mps 3           # split on variable 3
# =============================================================================

using CuttingPlane
using JuMP

include(joinpath(@__DIR__, "examples", "toy_model.jl"))

# ---- 1. Load the MIP --------------------------------------------------------
integer_model = isempty(ARGS) ? toy_model() : read_model(ARGS[1])
println(integer_model)

# ---- 2. Solve the MIP and its LP relaxation ---------------------------------
optimize_model!(integer_model)
z_int, _ = retrieve_objective_values(integer_model)

relaxed_model = relax_model(integer_model)
optimize_model!(relaxed_model)
z_lp, _ = retrieve_objective_values(relaxed_model)
ω, var_names = get_variable_values(relaxed_model)

print_banner("LP relaxation")
print_variable_values(ω, var_names)
println("LP bound: ", z_lp, "    MIP optimum: ", z_int)

integer_vars = integer_variable_indices(integer_model)
fractional = fractional_integer_variables(relaxed_model, integer_vars)
if isempty(fractional)
    println("\nThe LP optimum is already integral; no cut needed.")
    exit()
end

# ---- 3. Choose the split variable -------------------------------------------
var_index = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : first(fractional)
var_index in fractional ||
    error("Variable $var_index is not fractional. Fractional integer variables: $fractional")
println("Splitting on variable $var_index ($(var_names[var_index]) = $(ω[var_index]))")

# ---- 4. Build the two sides of the split disjunction ------------------------
standard_form_model = create_standard_form_model(relaxed_model)
subproblems, D0 = create_disjunctive_subproblems(standard_form_model, var_index, integer_vars)
for (t, sp) in enumerate(subproblems)
    optimize_model!(sp)
    if termination_status(sp) == MOI.INFEASIBLE
        # One side is empty, so the other side alone is the disjunctive hull:
        # adding its split row to the relaxation is already the best cut.
        println("\nSide $t of the disjunction is infeasible; the split row itself is the cut (no VPC needed).")
        exit()
    end
end

# ---- 5. Points and rays of each side, then the point-ray LP -----------------
results = [process_subproblem(sp) for sp in subproblems]
points = [r.point for r in results]
rays = [r.rays for r in results]
prlp, α, β = solve_prlp(ω, points, rays; normalization = :auto)
α' * ω < β - 1e-6 || println("Note: the cut does not separate the LP optimum (α'ω = $(α' * ω), β = $β).")

# ---- 6. Add the VPC and measure the improvement -----------------------------
print_banner("V-polyhedral cut")
vpc_model = add_cut(relaxed_model, α, β)
optimize_model!(vpc_model)
z_vpc, _ = retrieve_objective_values(vpc_model)
println("Bound with cut: ", z_vpc, "    gap closed: ", round(gap_closed(z_vpc, z_int, z_lp); digits = 2), "%")
print_variable_values(get_variable_values(vpc_model)...)

# ---- 7. Monoidal strengthening ----------------------------------------------
print_banner("Monoidal strengthened cut")
U = prepare_monoidal_data(subproblems, results, α)
α̃ = monoidal_strengthening(subproblems, var_index, U, α, integer_vars)
if α̃ ≈ α
    println("Strengthening leaves the cut unchanged: α = ", α)
else
    svpc_model = add_cut(relaxed_model, α̃, β)
    optimize_model!(svpc_model)
    z_svpc, _ = retrieve_objective_values(svpc_model)
    println("Bound with strengthened cut: ", z_svpc, "    gap closed: ",
            round(gap_closed(z_svpc, z_int, z_lp); digits = 2), "%")
    print_variable_values(get_variable_values(svpc_model)...)
end
