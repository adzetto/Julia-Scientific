# File: galerkin_fem_optimized.jl
# Description: Optimized solution for a boundary value problem using Galerkin FEM with linear basis functions
# Author: adzetto
# Date: 2024-06-20

#=
5. Use the Galerkin finite element method with continuous piecewise linear basis functions to solve the problem of Exercise 2, but with homogeneous Dirichlet boundary conditions:
$$
\begin{gathered}
\frac{d}{d x}\left(\left(1+x^2\right) \frac{d u}{d x}\right)=f(x), \quad 0 \leq x \leq 1, \\
u(0)=0, \quad u(1)=0 .
\end{gathered}
$$
(a) Derive the matrix equation that you will need to solve for this problem.
(b) Write a MATLAB code to solve this set of equations. You can test your code on a problem where you know the solution by choosing a function $u(x)$ that satisfies the boundary conditions and determining what $f(x)$ must be in order for $u(x)$ to solve the problem. Try $u(x)=x(1-x)$. Then $f(x)=-2\left(3 x^2-x+1\right)$.
(c) Try several different values for the mesh size $h$. Based on your results, what would you say is the order of accuracy of the Galerkin method with continuous piecewise linear basis functions?
=#

using LinearAlgebra
using SparseArrays
using Plots

# Stiffness matrix
function assemble_matrix(N)
    h = 1.0 / N
    diagonal = 2 / h .+ 2h / 3
    off_diagonal = -1 / h .+ h / 6
    
    A = spdiagm(-1 => fill(off_diagonal, N-2),
                0 => fill(diagonal, N-1),
                1 => fill(off_diagonal, N-2))
    
    return A
end

function assemble_rhs(N, f)
    h = 1.0 / N
    x = range(h, stop=1-h, length=N-1)
    b = h * (f.(x .- h/2) .+ f.(x .+ h/2)) / 2
    
    return b
end

# Galerkin FEM
function galerkin_fem(N, f)
    A = assemble_matrix(N)
    b = assemble_rhs(N, f)
    u = A \ b
    return [0; u; 0]  # B.C.
end

u_exact(x) = x * (1 - x)
f(x) = -2 * (3 * x^2 - x + 1)

function error_analysis(N_values)
    errors = []

    for N in N_values
        u_computed = galerkin_fem(N, f)
        x_values = range(0, stop=1, length=N+1)
        u_exact_values = u_exact.(x_values)
        error = norm(u_computed - u_exact_values, Inf)
        push!(errors, error)
    end

    for (N, error) in zip(N_values, errors)
        println("N = $N, Error = $error")
    end

    return N_values, errors
end

function main()
    N_values = [10, 20, 40, 80, 160]
    N_values, errors = error_analysis(N_values)
    
    p = plot(1.0 ./ N_values, errors, marker=:o, xscale=:log10, yscale=:log10,
    xlabel="Mesh size (h)", ylabel="Error", title="Error vs Mesh Size",
    label="Error", grid=true)
    savefig(p, "error_vs_mesh_size.pdf")
end

main()