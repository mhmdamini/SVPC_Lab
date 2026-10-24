"""
    CuttingPlane

V-polyhedral disjunctive cuts (VPCs) for mixed-integer linear programs, with an
optional monoidal strengthening step, built on JuMP and HiGHS.

Pipeline (see `run.jl` for an end-to-end script):

 1. Solve the MIP and its LP relaxation.
 2. Rewrite the relaxation in standard form (`max c'x, Ax <= b, x >= 0`).
 3. Pick a fractional integer variable `x_k` and build the split disjunction
    `x_k <= ⌊x̄_k⌋  ∨  x_k >= ⌈x̄_k⌉` (two subproblems).
 4. For each subproblem, extract the optimal vertex and the extreme rays of
    the LP cone at that vertex from the optimal simplex tableau.
 5. Solve the point-ray LP (PRLP) to obtain a valid cut `α'x >= β`.
 6. Add the cut to the relaxation, re-solve, and report how much of the
    integrality gap was closed.
 7. (Optional) Strengthen the cut coefficients with monoidal strengthening.

Gomory fractional cuts (`gomory_fractional_cut`), standard intersection cuts
(`sic_cut`) and Gomory mixed-integer cuts (`gmi_cut`) are also available, and `generate_cut`
produces any of these cut families from a solved relaxation.
"""
module CuttingPlane

using DataFrames
using HiGHS
using JuMP
using LinearAlgebra

const MOI = JuMP.MOI

include("utils.jl")     # reading, relaxing, solving and reporting
include("forms.jl")     # standard-form and normal-form (slack) reformulations
include("tableau.jl")   # simplex tableaux, basis information and cone rays
include("vpc.jl")       # disjunctions, PRLP and cut generation
include("monoidal.jl")  # monoidal strengthening of a VPC
include("simplex.jl")   # textbook two-phase simplex (iterations, cross-check)
include("tableau_cuts.jl") # Gomory fractional, intersection and GMI cuts
include("cuts.jl")      # generate_cut: one entry point for every cut family
include("iterative.jl") # iterative cutting-plane method (one or several families)
include("branch.jl")    # branch and bound / branch and cut

# utils.jl
export read_model, relax_model, integer_variable_indices, optimize_model!, get_variable_values,
       retrieve_objective_values, gap_closed, print_variable_values,
       print_banner
# forms.jl
export create_standard_form_model, create_normal_form_model
# tableau.jl
export generate_first_simplex_tableau, generate_optimal_simplex_tableau,
       compute_rays, verify_optimal_tableau
# simplex.jl
export tableau_simplex, SimplexStep
# vpc.jl
export find_fractional_elements, fractional_integer_variables,
       create_disjunctive_subproblems, process_subproblem, solve_prlp,
       add_cut
# monoidal.jl
export prepare_monoidal_data, monoidal_strengthening
# tableau_cuts.jl, cuts.jl
export gomory_fractional_cut, integer_scaled_model, sic_cut, gmi_cut,
       intersection_step_lengths, generate_cut, CUT_METHODS
# iterative.jl
export cutting_plane_loop, separation_round!, AddedCut, integer_cut
# branch.jl
export branch_and_bound, BBNode

end # module
