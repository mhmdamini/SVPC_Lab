# =============================================================================
# Textbook two-phase tableau simplex, independent of the LP solver.
#
# It solves the standard form  max c'x  s.t.  A x <= b, x >= 0  with slacks
# s (A x + s = b) and records every tableau, so the iterations can be shown
# step by step and the HiGHS-based optimal tableau can be cross-checked.
#
# * Phase I  (only if some b_i < 0): row i is negated (−A_i x − s_i = −b_i)
#   and gets an artificial a_i; maximise −Σ a. A negative optimum means the
#   LP is infeasible. Remaining basic artificials (at zero) are pivoted out,
#   or their row is dropped as redundant.
# * Phase II: maximise c'x from the feasible basis.
# * Bland's rule (smallest index enters; ties in the ratio test leave by
#   smallest index) guarantees termination, also on degenerate problems.
#
# Tableau convention, as in tableau.jl: rows  x_B + B⁻¹N x_N = B⁻¹b  and the
# objective row  z + d'x = z̄  with reduced costs d = c_B'B⁻¹A − c.
# =============================================================================

"""One recorded simplex step (the tableau *before* the pivot is applied)."""
struct SimplexStep
    phase::Int                       # 1 or 2
    columns::Vector{String}          # variable names, in column order
    basis::Vector{String}            # basic variable of each row
    body::Matrix{Float64}            # constraint rows  (m × ncols)
    rhs::Vector{Float64}             # B⁻¹b
    reduced::Vector{Float64}         # objective row d
    objective::Float64               # z̄ (phase objective)
    entering::Union{Nothing,Int}     # pivot column, or nothing (optimal / final)
    leaving::Union{Nothing,Int}      # pivot row
    ratios::Vector{Union{Nothing,Float64}}   # ratio test per row (nothing = not eligible)
    note::String
end

"""
    tableau_simplex(std_model; max_iter = 500, tol = 1e-9)

Run the two-phase tableau simplex on a standard-form model (as produced by
`create_standard_form_model`). Returns a named tuple

* `status`     – `:optimal`, `:infeasible`, `:unbounded` or `:iteration_limit`;
* `steps`      – `Vector{SimplexStep}`; the last step is the final tableau;
* `x`          – structural solution (when optimal), `objective` – c'x + offset;
* `basis`      – final basic variables; `columns` – names of all columns.
"""
function tableau_simplex(std_model::Model; max_iter::Int = 500, tol::Real = 1e-9)
    data = lp_matrix_data(std_model)
    A, b, c = Matrix(data.A), copy(data.b_upper), data.c
    m, n = size(A)
    xnames = [isempty(name(v)) ? "x$j" : name(v) for (j, v) in enumerate(data.variables)]
    snames = ["s$i" for i in 1:m]

    # Rows with b_i < 0 are negated and receive an artificial variable.
    neg = findall(<(0), b)
    k = length(neg)
    T = hcat(A, Matrix(1.0I, m, m), zeros(m, k))
    for (a, i) in enumerate(neg)
        T[i, :] .*= -1
        b[i] = -b[i]
        T[i, n+m+a] = 1.0
    end
    columns = vcat(xnames, snames, ["a$i" for i in neg])
    basis = [i in neg ? n + m + findfirst(==(i), neg) : n + i for i in 1:m]   # column index per row
    steps = SimplexStep[]

    # Objective row  d = c_B'·T − c_full,  value z̄ = c_B'·b.
    objective_row(cfull) = (cfull[basis]' * T)' .- cfull, cfull[basis]' * b

    function record!(phase, cfull, entering, leaving, ratios, note)
        d, z = objective_row(cfull)
        push!(steps, SimplexStep(phase, copy(columns), columns[basis], copy(T), copy(b), d,
                                 z + (phase == 2 ? data.c_offset : 0.0), entering, leaving, ratios, note))
    end

    # Run simplex iterations on the current tableau for objective `cfull`.
    function run_phase!(phase, cfull)
        for _ in 1:max_iter
            d, _ = objective_row(cfull)
            j = findfirst(<(-tol), d)                         # Bland: smallest improving index
            if j === nothing
                record!(phase, cfull, nothing, nothing, fill(nothing, length(b)),
                        phase == 1 ? "Phase I optimal." : "Optimal: no reduced cost is negative.")
                return :optimal
            end
            ratios = Union{Nothing,Float64}[T[i, j] > tol ? b[i] / T[i, j] : nothing for i in eachindex(b)]
            eligible = findall(!isnothing, ratios)
            if isempty(eligible)
                record!(phase, cfull, j, nothing, ratios, "Unbounded: column $(columns[j]) has no positive entry.")
                return :unbounded
            end
            best = minimum(ratios[i] for i in eligible)
            ties = filter(i -> ratios[i] <= best + tol, eligible)
            r = ties[argmin(basis[ties])]                     # Bland: smallest basic index leaves
            record!(phase, cfull, j, r, ratios,
                    "$(columns[j]) enters (reduced cost $(round(d[j]; sigdigits = 4))), " *
                    "$(columns[basis[r]]) leaves (min ratio $(round(best; sigdigits = 4))).")
            # Pivot on T[r, j].
            p = T[r, j]
            T[r, :] ./= p; b[r] /= p
            for i in axes(T, 1)
                i == r && continue
                f = T[i, j]
                abs(f) > 0 || continue
                T[i, :] .-= f .* T[r, :]; b[i] -= f * b[r]
            end
            basis[r] = j
        end
        return :iteration_limit
    end

    # Value of structural variable j in the current basic solution.
    structural_value(j) = (r = findfirst(==(j), basis); r === nothing ? 0.0 : b[r])
    function result(status)
        x = status == :optimal ? [structural_value(j) for j in 1:n] : nothing
        return (status = status, steps = steps, columns = columns, basis = columns[basis], x = x,
                objective = x === nothing ? nothing : c' * x + data.c_offset)
    end

    # ---- Phase I ----
    if k > 0
        c1 = vcat(zeros(n + m), -ones(k))
        st = run_phase!(1, c1)
        st == :optimal || return result(st)
        if steps[end].objective < -1e-7
            steps[end] = SimplexStep(1, steps[end].columns, steps[end].basis, steps[end].body, steps[end].rhs,
                                     steps[end].reduced, steps[end].objective, nothing, nothing, steps[end].ratios,
                                     "Infeasible: Phase I optimum $(round(steps[end].objective; sigdigits = 4)) < 0.")
            return result(:infeasible)
        end
        # Drive remaining (zero-valued) artificials out of the basis.
        for r in reverse(eachindex(basis))
            basis[r] > n + m || continue
            j = findfirst(jj -> abs(T[r, jj]) > tol, 1:n+m)
            if j === nothing                                   # redundant row
                T = T[setdiff(axes(T, 1), r), :]; b = b[setdiff(eachindex(b), r)]; deleteat!(basis, r)
            else
                p = T[r, j]; T[r, :] ./= p; b[r] /= p
                for i in axes(T, 1)
                    i == r && continue
                    f = T[i, j]; T[i, :] .-= f .* T[r, :]; b[i] -= f * b[r]
                end
                basis[r] = j
            end
        end
        T = T[:, 1:n+m]                                        # drop artificial columns
        columns = columns[1:n+m]
    end

    # ---- Phase II ----
    st = run_phase!(2, vcat(c, zeros(m)))
    return result(st)
end
