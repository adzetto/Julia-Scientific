# linear_interpolation_methods.jl
# Description:
# Applys linear interpolation and the secant technique to solve a number of equations.
# Author: adzetto
# Date: 2024-06-18

using Printf

function linear_interpolation(f, x0, x1, tol, max_iter)
    for iter in 1:max_iter
        fx0 = f(x0)
        fx1 = f(x1)
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        fx2 = f(x2)
        if abs(fx2) < tol
            return x2, iter
        end
        x0, x1 = x1, x2
    end
    return x1, max_iter
end

function secant_method(f, x0, x1, tol, max_iter)
    for iter in 1:max_iter
        fx0 = f(x0)
        fx1 = f(x1)
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        if abs(fx1 - fx0) < tol
            return x2, iter
        end
        x0, x1 = x1, x2
    end
    return x1, max_iter
end

# Exercise 5: Polynomial root approximation
f1(x) = x^3 + x^2 - 3x - 3
root1, iter1 = linear_interpolation(f1, -2, -1, 1e-6, 1000)
@printf("Root of x^3 + x^2 - 3x - 3 using linear interpolation: %.8f (Iterations: %d)\n", root1, iter1)

# Exercise 6: Secant method for different starting values
root2a, iter2a = secant_method(f1, -1.5, -1.7, 1e-6, 1000)
root2b, iter2b = secant_method(f1, -1.5, -1.1, 1e-6, 1000)
root2c, iter2c = secant_method(f1, -1.5, -1.25, 1e-6, 1000)
@printf("Root using secant method with start values -1.5, -1.7: %.8f (Iterations: %d)\n", root2a, iter2a)
@printf("Root using secant method with start values -1.5, -1.1: %.8f (Iterations: %d)\n", root2b, iter2b)
@printf("Root using secant method with start values -1.5, -1.25: %.8f (Iterations: %d)\n", root2c, iter2c)

# Exercise 7: Intersection of y = x^3 - x + 1 and y = 2x^2
f2(x) = x^3 - x + 1 - 2x^2
root3, iter3 = secant_method(f2, -2, 2, 1e-6, 1000)
@printf("Intersection points of y = x^3 - x + 1 and y = 2x^2: %.8f (Iterations: %d)\n", root3, iter3)

# Exercise 8a: Linear interpolation for provided equations
f3(x) = x^2 - 4
root4, iter4 = linear_interpolation(f3, 0, 3, 1e-6, 1000)
@printf("Root of f3 using linear interpolation: %.8f (Iterations: %d)\n", root4, iter4)

# Exercise 8b: Modified linear interpolation
# One possible modification is to enhance convergence by using averaging.
function modified_linear_interpolation(f, x0, x1, tol, max_iter)
    for iter in 1:max_iter
        fx0 = f(x0)
        fx1 = f(x1)
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        fx2 = f(x2)
        if abs(fx2) < tol
            return x2, iter
        end
        x0, x1 = x1, (x1 + x2) / 2  # Modify: Calculate the average of x1 and x2 to improve convergence.
    end
    return x1, max_iter
end

root5, iter5 = modified_linear_interpolation(f3, 0, 3, 1e-6, 1000)
@printf("Root of f3 using modified linear interpolation: %.8f (Iterations: %d)\n", root5, iter5)

# Exercise 9: Algorithm for the secant method (defined above as secant_method)