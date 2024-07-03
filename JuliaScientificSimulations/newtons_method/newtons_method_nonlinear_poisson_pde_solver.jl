# newtons_method_nonlinear_poisson_pde_solver.jl
# Description:
# This script solves a nonlinear Poisson PDE using Newton's method.
# Author: adzetto
# Date: 18/06/2024


#=
4. Consider the nonlinear PDE in two dimensions given by
$$
-\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)+e^u=g(x, y)
$$

Here $u=u(x, y)$ is a function in the two variables $x$ and $y$, defined on the unit square $(0,1) \times(0,1)$ and subject to homogeneous Dirichlet boundary conditions. We discretize it on a uniform mesh, just like in Example 7.1, using $h=1 /(N+1)$, and obtain the equations
$$
\begin{aligned}
4 u_{i, j}-u_{i+1, j}-u_{i-1, j}-u_{i, j+1}-u_{i, j-1}+h^2 e^{u_{i, j}} & =h^2 g_{i, j}, \quad 1 \leq i, j \leq N \\
u_{i, j} & =0 \text { otherwise. }
\end{aligned}
$$
(a) Find the Jacobian matrix, $J$, and show that it is always symmetric positive definite.
(b) Write a JULIA program that solves this system of nonlinear equations using Newton's method for $N=8,16,32$. Generate a right-hand-side function $g(x, y)$ such that the exact solution of the differential problem is $u(x, y) \equiv 1$. Start with an initial guess of all zeros, and stop the iteration when $\left\|\delta u^{(k)}\right\|_2<10^{-6}$. Plot the norms $\left\|\delta u^{(k)}\right\|_2$ and explain the convergence behavior you observe.
=#

using LinearAlgebra
using SparseArrays
using Plots

N_values = [8, 16, 32]
tol = 1e-6
max_iter = 100

u_exact(x, y) = 1
g(x, y) = exp(1) + 2*π^2

function h(N)
    return 1 / (N + 1)
end

function F(u, N, h)
    F = zeros(N, N)
    for i in 1:N, j in 1:N
        F[i, j] = 4 * u[i, j] - (i < N ? u[i+1, j] : 0) - (i > 1 ? u[i-1, j] : 0) - (j < N ? u[i, j+1] : 0) - (j > 1 ? u[i, j-1] : 0) + h^2 * exp(u[i, j]) - h^2 * g(i*h, j*h)
    end
    return F
end

function J(u, N, h)
    J = spzeros(N*N, N*N)
    for i in 1:N, j in 1:N
        index = (i-1)*N + j
        J[index, index] = 4 + h^2 * exp(u[i, j])
        if i < N
            J[index, index+N] = -1
        end
        if i > 1
            J[index, index-N] = -1
        end
        if j < N
            J[index, index+1] = -1
        end
        if j > 1
            J[index, index-1] = -1
        end
    end
    return J
end

function newtons_method(F, J, N, h, tol, max_iter)
    u = zeros(N, N)
    norm_deltas = []
    
    for k in 1:max_iter
        F_val = F(u, N, h)
        J_val = J(u, N, h)
        delta_u = -J_val \ reshape(F_val, N*N)
        delta_u = reshape(delta_u, N, N)
        u += delta_u
        norm_delta = norm(delta_u, 2)
        push!(norm_deltas, norm_delta)
        
        if norm_delta < tol
            println("Converged in $k iterations.")
            return u, norm_deltas
        end
    end
    
    println("Did not converge within the maximum number of iterations.")
    return u, norm_deltas
end

for N in N_values
    h_val = h(N)
    u, norm_deltas = newtons_method(F, J, N, h_val, tol, max_iter)
    
    plot(norm_deltas, yscale=:log10, xlabel="Iteration", ylabel="||Δu^(k)||_2", title="Newton's Method Convergence for N=$N")
    savefig("/mnt/c/Users/lenovo/Documents/GitHub/Julia-Scientific/newton_convergence_N_$N.pdf")
end
