# minimization_newton_bfgs.jl
# Description:
# This script minimizes the function φ(x) = c^T x + 1/2 x^T H x using both Newton and BFGS methods.
# Author: adzetto
# Date: 18/06/2024

#=
13. Consider minimizing the function $\phi(\mathbf{x})=\mathbf{c}^T \mathbf{x}+\frac{1}{2} \mathbf{x}^T H \mathbf{x}$, where $\mathbf{c}=(5.04,-59.4,146.4,-96.6)^T$ and
$$
H=\left(\begin{array}{cccc}
.16 & -1.2 & 2.4 & -1.4 \\
-1.2 & 12.0 & -27.0 & 16.8 \\
2.4 & -27.0 & 64.8 & -42.0 \\
-1.4 & 16.8 & -42.0 & 28.0
\end{array}\right) .
$$

Try both Newton and BFGS methods, starting from $\mathbf{x}_0=(-1,3,3,0)^T$. Explain why the BFGS method requires significantly more iterations than Newton's.
=#

using LinearAlgebra
using Optim

c = [5.04, -59.4, 146.4, -96.6]
H = [0.16 -1.2 2.4 -1.4;
     -1.2 12.0 -27.0 16.8;
     2.4 -27.0 64.8 -42.0;
     -1.4 16.8 -42.0 28.0]

function φ(x)
    return dot(c, x) + 0.5 * dot(x, H * x)
end

function ∇φ(x)
    return c + H * x
end

function ∇∇φ(x)
    return H
end

x0 = [-1.0, 3.0, 3.0, 0.0]

newton_result = optimize(φ, ∇φ, ∇∇φ, x0, Newton(); inplace=false)
println("Newton's method result:")
println("Minimum value: ", newton_result.minimum)
println("Minimizer: ", newton_result.minimizer)
println("Iterations: ", newton_result.iterations)

bfgs_result = optimize(φ, ∇φ, x0, BFGS(); inplace=false)
println("\nBFGS method result:")
println("Minimum value: ", bfgs_result.minimum)
println("Minimizer: ", bfgs_result.minimizer)
println("Iterations: ", bfgs_result.iterations)
