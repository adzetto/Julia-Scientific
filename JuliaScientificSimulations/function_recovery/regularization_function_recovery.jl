# regularization_function_recovery.jl
# Description:
# This script recovers a function from noisy data using regularization.
# Author: adzetto
# Date: 18/06/2024

#=
17. This exercise is concerned with recovering a function $u(t)$ on the interval $[0,1]$ given noisy data $b_i$ at points $t_i=i h, i=0,1, \ldots, N$, with $N=1 / h$. Because the data values are noisy, we cannot simply set $u\left(t_i\right) \equiv u_i=b_i$ : knowing that $u(t)$ should be piecewise smooth, we add a regularization term to penalize excessive roughness in $u$. For the unknown vector $\mathbf{u}=$ $\left(u_0, u_1, \ldots, u_N\right)^T$ we therefore solve
$$
\min \phi_2(\mathbf{u})=\frac{h}{2} \sum_{i=1}^N \frac{1}{2}\left[\left(u_i-b_i\right)^2+\left(u_{i-1}-b_{i-1}\right)^2\right]+\frac{\beta h}{2} \sum_{i=1}^N\left(\frac{u_i-u_{i-1}}{h}\right)^2 .
$$
(a) Write down the gradient and the Hessian of this objective function. To describe the regularization term use the matrix
$$
W=\frac{1}{\sqrt{h}}\left(\begin{array}{ccccc}
-1 & 1 & & & \\
& -1 & 1 & & \\
& & \ddots & \ddots & \\
& & & -1 & 1
\end{array}\right) \in \Re^{N \times(N+1)}
$$

(b) Solve this problem numerically for the following problem instances. To "synthesize data" for a given $N$, start with
$$
b_p(t)= \begin{cases}1, & 0 \leq t<.25, \\ 2, & .25 \leq t<.5, \\ 2-100(t-.5)(.7-t), & .5 \leq t<.7, \\ 4, & .7 \leq t \leq 1\end{cases}
$$

Evaluate this at the grid points and then add noise as follows:
$$
\begin{aligned}
& \text { noisev = randn(size(b_p)) * mean(abs (b_p ) ) * noise; } \\
& \text { data = b_p + noisev; }
\end{aligned}
$$

The resulting values are the data that your program "sees." An illustration is given in Figure 9.13.

Figure 9.13. A depiction of the noisy function (in solid blue) and the function to be recovered (in red) for Exercise 17.
Plot this data and the recovered curve $\mathbf{u}$ for the parameter values $(\beta$, noise $)=\left(10^{-3}, .01\right)$, $\left(10^{-3}, .1\right),\left(10^{-4}, .01\right)$ and $\left(10^{-4}, .1\right)$. Try $N=64$ or $N=128$. What are your observations?

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

function φ2(u, b, h, β, W)
    return (h / 2) * sum((u .- b).^2) + (β / 2) * norm(W * u)^2
end

function ∇φ2(u, b, h, β, W)
    return h * (u .- b) + β * (W' * W) * u
end

function ∇∇φ2(h, β, W)
    return h * I + β * (W' * W)
end

W = spzeros(N, N+1)
for i in 1:N
    W[i, i] = -1 / sqrt(h)
    W[i, i+1] = 1 / sqrt(h)
end

function solve_regularization(b, h, β, W, tol=1e-6, max_iter=1000)
    n = length(b)
    u = zeros(n)
    H = ∇∇φ2(h, β, W)
    
    for k in 1:max_iter
        grad = ∇φ2(u, b, h, β, W)
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

parameter_values = [(1e-3, 0.01), (1e-3, 0.1), (1e-4, 0.01), (1e-4, 0.1)]
b_values = add_noise(b_p_values, 0.1)

plot(t, b_p_values, label="Original function", linewidth=2)
plot!(t, b_values, label="Noisy data", linestyle=:dash)

for (β, noise_level) in parameter_values
    b_noisy = add_noise(b_p_values, noise_level)
    u_recovered = solve_regularization(b_noisy, h, β, W)
    plot!(t, u_recovered, label="Recovered (β=$β, noise=$noise_level)")
end

xlabel!("t")
ylabel!("u(t)")
title!("Function Recovery with Regularization")
savefig("/mnt/c/Users/lenovo/Documents/GitHub/Julia-Scientific/function_recovery.pdf")
