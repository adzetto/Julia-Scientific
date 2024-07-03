# regularization_function_recovery_iterated_least_squares.jl
# Description:
# This script recovers a function from noisy data using iterated least squares with regularization.
# Author: adzetto
# Date: 18/06/2024

#=
18. Continuing Exercise 17, unfortunately if the data function contains jump discontinuities, then these discontinuities will be smeared by the regularization proposed there. So, consider instead the minimization problem
$$
\min \phi_1(\mathbf{u})=\frac{h}{2} \sum_{i=1}^N \frac{1}{2}\left[\left(u_i-b_i\right)^2+\left(u_{i-1}-b_{i-1}\right)^2\right]+\gamma h \sum_{i=1}^N \sqrt{\left(\frac{u_i-u_{i-1}}{h}\right)^2+\varepsilon}
$$
where $\varepsilon=10^{-6}$, say.

(a) Let $T_{1, h}$ denote what $\gamma$ multiplies in the objective function $\phi_1$. Show that
$$
\frac{\partial T_{1, h}}{\partial u_j}=\frac{u_j-u_{j-1}}{\sqrt{\left(u_j-u_{j-1}\right)^2+\varepsilon h^2}}+\frac{u_j-u_{j+1}}{\sqrt{\left(u_j-u_{j+1}\right)^2+\varepsilon h^2}} .
$$

Moreover, letting
$$
\hat{D}=\operatorname{diag}\left\{h / \sqrt{\left(u_j-u_{j-1}\right)^2+\varepsilon h^2}\right\}, \quad \hat{B}=\sqrt{\hat{D}},
$$
it is possible to write
$$
\nabla T_{1, h}=W^T \hat{B}^T \hat{B} W \mathbf{u}=W^T \hat{D} W \mathbf{u} .
$$
(b) A method for solving this problem then suggests itself, whereby at the beginning of each iteration we fix $\hat{D}$ based on the current iterate and apply the usual algorithms for a linear weighted least squares functional. This is called iterated least squares.
Use this method to solve the same problems as in the previous question (i.e., same synthesized data), for $\gamma=10^{-2}$ and $\gamma=10^{-1}$. Discuss your observations.
=#

using LinearAlgebra
using SparseArrays
using Random
using Plots
using Statistics

function b_p(t)
    if 0 <= t < 0.25
        return 1.0
    elseif 0.25 <= t < 0.5
        return 2.0
    elseif 0.5 <= t < 0.7
        return 2.0 - 100.0 * (t - 0.5) * (0.7 - t)
    else
        return 4.0
    end
end

N = 128
h = 1.0 / (N + 1)
t = [i * h for i in 0:N]
b_p_values = [b_p(ti) for ti in t]

function add_noise(data, noise_level)
    noise = randn(length(data)) * mean(abs.(data)) * noise_level
    return data .+ noise
end

function φ1(u, b, h, γ, W, ϵ)
    return (h / 2) * sum((u .- b).^2) + γ * h * sum(sqrt.((W * u).^2 .+ ϵ))
end

function ∇φ1(u, b, h, γ, W, ϵ)
    D_inv = diagm(1.0 ./ sqrt.((W * u).^2 .+ ϵ))
    return h * (u .- b) + γ * W' * D_inv * W * u
end

function ∇∇φ1(u, h, γ, W, ϵ)
    D_inv = diagm(1.0 ./ sqrt.((W * u).^2 .+ ϵ))
    return h * I + γ * W' * D_inv * W
end

W = spzeros(N, N+1)
for i in 1:N
    W[i, i] = -1 / sqrt(h)
    W[i, i+1] = 1 / sqrt(h)
end

function solve_iterated_least_squares(b, h, γ, W, ϵ, tol=1e-6, max_iter=1000)
    n = length(b)
    u = zeros(n)
    
    for k in 1:max_iter
        H = ∇∇φ1(u, h, γ, W, ϵ)
        grad = ∇φ1(u, b, h, γ, W, ϵ)
        delta_u = -H \ grad
        u += delta_u
        
        if norm(delta_u) < tol
            println("Converged in $k iterations.")
            return u
        end
    end
    
    println("Did not converge within the maximum number of iterations.")
    return u
end

parameter_values = [(1e-2, 0.01), (1e-2, 0.1), (1e-1, 0.01), (1e-1, 0.1)]
ϵ = 1e-6
b_values = add_noise(b_p_values, 0.1)

plot(t, b_p_values, label="Original function", linewidth=2)
plot!(t, b_values, label="Noisy data", linestyle=:dash)

for (γ, noise_level) in parameter_values
    b_noisy = add_noise(b_p_values, noise_level)
    u_recovered = solve_iterated_least_squares(b_noisy, h, γ, W, ϵ)
    plot!(t, u_recovered, label="Recovered (γ=$γ, noise=$noise_level)")
end

xlabel!("t")
ylabel!("u(t)")
title!("Function Recovery with Iterated Least Squares")
savefig("/mnt/c/Users/lenovo/Documents/GitHub/Julia-Scientific/function_recovery_iterated_least_squares.pdf")
