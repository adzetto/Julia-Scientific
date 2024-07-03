# newton_iteration_a.jl
# Description:
# In order to resolve the set of nonlinear equations,
#this script applies the Newton-Raphson algorithm.
#              x1^2 + x1*x2^3 = 9
#              3*x1^2*x2 - x2^3 = 4
# Author: adzetto
# Date: 18/06/2024


using LinearAlgebra

function newton_iteration_a(x_init, tol=1e-8, max_iter=100)
    function F(x)
        [x[1]^2 + x[1] * x[2]^3 - 9;
         3 * x[1]^2 * x[2] - x[2]^3 - 4]
    end

    function J(x)
        [2 * x[1] + x[2]^3  3 * x[1] * x[2]^2;
         6 * x[1] * x[2]  3 * x[1]^2 - 3 * x[2]^2]
    end

    x = x_init
    for k in 1:max_iter
        Fx = F(x)
        Jx = J(x)
        delta = Jx \ Fx
        x = x - delta
        if norm(delta) < tol
            break
        end
    end
    return x
end

x_init = [1.0, 1.0]
solution, iterations = newton_iteration_a(x_init)
println("Solution: ", solution)
println("Iterations: ", iterations)

# newton_iteration_b.jl
# Description:
# In order to resolve the set of nonlinear equations,
#this script applies the Newton-Raphson algorithm.
#              x1 + x2 - 2*x1*x2 = 0
#              x1^2 + x2^2 - 2*x1 + 2*x2 = -1
# Author: adzetto
# Date: 18/06/2024

using LinearAlgebra

function newton_iteration_b(x_init, tol=1e-8, max_iter=100)
    function F(x)
        return [x[1] + x[2] - 2 * x[1] * x[2];
                x[1]^2 + x[2]^2 - 2 * x[1] + 2 * x[2] + 1]
    end

    function J(x)
        return [1 - 2 * x[2]  1 - 2 * x[1];
                2 * x[1] - 2  2 * x[2] + 2]
    end

    x = x_init
    for k in 1:max_iter
        Fx = F(x)
        Jx = J(x)
        delta = Jx \ Fx
        x = x - delta
        if norm(delta) < tol
            return x, k
        end
    end
    error("Newton iteration did not converge")
end

x_init = [1.0, 1.0]
solution, iterations = newton_iteration_b(x_init)
println("Solution: ", solution)
println("Iterations: ", iterations)


# newton_iteration_c.jl
# Description:
# In order to resolve the set of nonlinear equations,
#this script applies the Newton-Raphson algorithm.
#              x1^3 - x2^2 = 0
#              x1 + x1^2*x2 = 2
# Author: adzetto
# Date: 18/06/2024

using LinearAlgebra

function newton_iteration_c(x_init, tol=1e-8, max_iter=100)
    function F(x)
        return [x[1]^3 - x[2]^2;
                x[1] + x[1]^2 * x[2] - 2]
    end

    function J(x)
        return [3 * x[1]^2  -2 * x[2];
                1 + 2 * x[1] * x[2]  x[1]^2]
    end

    x = x_init
    for k in 1:max_iter
        Fx = F(x)
        Jx = J(x)
        delta = Jx \ Fx
        x = x - delta
        if norm(delta) < tol
            return x, k
        end
    end
    error("Newton iteration did not converge")
end

x_init = [1.0, 1.0]
solution, iterations = newton_iteration_c(x_init)
println("Solution: ", solution)
println("Iterations: ", iterations)
