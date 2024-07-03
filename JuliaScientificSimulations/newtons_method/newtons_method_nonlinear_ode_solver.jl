# newtons_method_nonlinear_ode_solver.jl
# Description:
# This script solves a nonlinear ODE using Newton's method.
# Author: adzetto
# Date: 18/06/2024

#=
7. Use Newton's method to solve a discretized version of the differential equation
$$
y^{\prime \prime}=-\left(y^{\prime}\right)^2-y+\ln x, \quad 1 \leq x \leq 2, \quad y(1)=0, \quad y(2)=\ln 2 .
$$

The discretization on a uniform mesh, with the notation of Example 9.3, can be
$$
\frac{y_{i+1}-2 y_i+y_{i-1}}{h^2}+\left(\frac{y_{i+1}-y_{i-1}}{2 h}\right)^2+y_i=\ln (i h), \quad i=1,2, \ldots, n
$$

The actual solution of this problem is $y(x)=\ln x$. Compare your numerical results to the solution $y(x)$ for $n=8,16,32$, and 64 . Make observations regarding the convergence behavior of Newton's method in terms of the iterations and the mesh size, as well as the solution error.
=#

using LinearAlgebra
using SparseArrays
using Plots

y_exact(x) = log(x)

# Define the discretized ODE function
function F(y, h, n)
    F = zeros(n)
    for i in 2:n-1
        F[i] = (y[i+1] - 2*y[i] + y[i-1]) / h^2 + ((y[i+1] - y[i-1]) / (2*h))^2 + y[i] - log(1 + i * h)
    end
    return F
end

# Define the Jacobian of the discretized ODE
function J(y, h, n)
    J = spzeros(n, n)
    for i in 2:n-1
        J[i, i-1] = 1/h^2 - (y[i+1] - y[i-1]) / (4*h^2)
        J[i, i] = -2/h^2 + 1 + (y[i+1] - y[i-1])^2 / (4*h^2)
        J[i, i+1] = 1/h^2 + (y[i+1] - y[i-1]) / (4*h^2)
    end
    return J
end

# Newton's method implementation
function newtons_method(F, J, y, h, n, tol, max_iter)
    norm_deltas = []
    lambda = 1e-6  # Regularization parameter

    for k in 1:max_iter
        F_val = F(y, h, n)
        J_val = J(y, h, n)
        # Add a small value to the diagonal of J to prevent it from being singular
        J_reg = J_val + lambda * I
        delta_y = -J_reg \ F_val
        y += delta_y
        norm_delta = norm(delta_y, 2)
        push!(norm_deltas, norm_delta)

        if norm_delta < tol
            println("Converged in $k iterations.")
            return y, norm_deltas
        end
    end

    println("Did not converge within the maximum number of iterations.")
    return y, norm_deltas
end

n_values = [8, 16, 32, 64]
tol = 1e-6
max_iter = 100
plot_combined = plot(title="Newton's Method Convergence", xlabel="Iteration", ylabel="||Δy^(k)||_2", yscale=:log10)

for n in n_values
    h = (2 - 1) / (n + 1)
    y = zeros(n+2)
    y[1] = 0  # Boundary condition y(1) = 0
    y[end] = log(2)  # Boundary condition y(2) = log(2)
    
    y, norm_deltas = newtons_method(F, J, y, h, n+2, tol, max_iter)
    
    plot!(plot_combined, norm_deltas, label="n=$(n+2)")
    
    x_values = [1 + i * h for i in 0:n+1]
    y_exact_values = y_exact.(x_values)
    println("For n = $(n+2):")
    println("Numerical solution: $y")
    println("Exact solution: $y_exact_values")
    println("Error: $(norm(y - y_exact_values, 2))")
end

savefig(plot_combined, "/mnt/c/Users/lenovo/Documents/GitHub/Julia-Scientific/newton_convergence_combined.pdf")
