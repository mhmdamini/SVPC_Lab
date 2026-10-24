# =============================================================================
# Monoidal strengthening of a V-polyhedral cut (Balas–Jeroslow, 1980) for the
# split disjunction  x_k <= π0  ∨  x_k >= π0 + 1.
#
# Side t is the standard-form system  A_t x <= b_t, x >= 0,  whose last row is
# the split row  d_t'x <= D_t  (d_1 = e_k, d_2 = -e_k). At the optimal vertex
# p_t of side t the VPC  α'x >= β  is valid for the basis cone, so by Farkas
#
#     α = -A_t'u_t + v_t,      α'p_t = -b_t'u_t >= β,      u_t, v_t >= 0,
#
# where u_t[i] = α'r for the ray r spanned by the slack of row i (0 if that
# slack is basic) and v_t[j] = α'r for the ray of a non-basic structural x_j.
# Dropping v_t gives the side-t coefficients  α_t = -A_t'u_t <= α.
#
# For every integer vector μ supported on the integer variables, π = e_k + μ
# defines another valid split  π'x <= π0 ∨ π'x >= π0 + 1, and replacing the
# split rows by these only changes the coefficients:
#
#     side 1:  α_1 - u0_1 μ,       side 2:  α_2 + u0_2 μ        (u0_t = u_t[end])
#
# Both sides keep the right-hand side -b_t'u_t >= β, so  α̃'x >= β  is valid
# whenever α̃ >= α_t(μ) on both sides. Coordinate-wise this gives, for integer j,
#
#     α̃_j = min_{μ ∈ ℤ} max(α_1j - u0_1 μ, α_2j + u0_2 μ)      (<= α_j at μ = 0)
#
# and α̃_j = α_j for continuous j.
# =============================================================================

"""
    prepare_monoidal_data(subproblems, results, α) -> U

Farkas multipliers of the rows of each side of the split for the cut `α'x >= β`.
For side `t`, `U[t][i]` is `α'r` for the ray `r` spanned by the slack of row
`i` when that slack is non-basic, and `0` otherwise (negative round-off is
clipped to 0). Rows follow `lp_matrix_data(subproblems[t])`, so the last entry
is the multiplier of the split row. `results[t]` is the output of
`process_subproblem(subproblems[t])`.
"""
function prepare_monoidal_data(subproblems, results, α::AbstractVector)
    U = Vector{Vector{Float64}}(undef, length(subproblems))
    for (t, res) in enumerate(results)
        n_rows = size(lp_matrix_data(subproblems[t]).A, 1)
        U[t] = zeros(n_rows)
        weights = α' * res.rays                       # α'r for every ray of side t
        n_struct = size(res.rays, 1)
        slack_names = name.(all_variables(res.normal_model))[n_struct+1:end]
        for (j, v) in enumerate(res.nonbasic)
            i = findfirst(==(v), slack_names)          # row whose slack spans ray j
            i === nothing || (U[t][i] = max(weights[j], 0.0))
        end
    end
    return U
end

"""
    best_shift(a1, u1, a2, u2; atol = 1e-9) -> Float64

`min_{μ ∈ ℤ} max(a1 - u1 μ, a2 + u2 μ)` for `u1, u2 >= 0`. The function is
convex and piecewise linear in μ, so the minimum is at the floor or ceiling of
the crossing point `(a1 - a2) / (u1 + u2)`.
"""
function best_shift(a1, u1, a2, u2; atol::Real = 1e-9)
    u1 = u1 > atol ? u1 : 0.0
    u2 = u2 > atol ? u2 : 0.0
    u1 + u2 == 0 && return max(a1, a2)
    μ = (a1 - a2) / (u1 + u2)
    f(m) = max(a1 - u1 * m, a2 + u2 * m)
    return min(f(floor(μ)), f(ceil(μ)))
end

"""
    monoidal_strengthening(subproblems, var_index, U, α, integer_vars) -> α̃

Monoidally strengthened coefficients `α̃ <= α` of the VPC `α'x >= β` for the
split on `x_{var_index}` (the right-hand side β is unchanged). Only the
coefficients of `integer_vars` are changed.

* `subproblems` – the two standard-form sides from `create_disjunctive_subproblems`
  (side 1 is `x_k <= π0`, side 2 is `x_k >= π0 + 1`), split row last;
* `U` – row multipliers from `prepare_monoidal_data`.
"""
function monoidal_strengthening(subproblems, var_index::Int, U, α::AbstractVector, integer_vars)
    length(subproblems) == 2 || error("Monoidal strengthening is implemented for split disjunctions only.")
    n = length(α)
    αt = Vector{Vector{Float64}}(undef, 2)
    for t in 1:2
        A = Matrix(lp_matrix_data(subproblems[t]).A)
        split_row = zeros(n)
        split_row[var_index] = t == 1 ? 1.0 : -1.0
        A[end, :] ≈ split_row || error("The last row of side $t is not the split row on variable $var_index.")
        αt[t] = -A' * U[t]                             # α = α_t + v_t with v_t >= 0
    end
    u0 = [last(U[1]), last(U[2])]

    α̃ = float.(collect(α))
    for j in integer_vars
        α̃[j] = min(α[j], best_shift(αt[1][j], u0[1], αt[2][j], u0[2]))
    end
    return α̃
end
