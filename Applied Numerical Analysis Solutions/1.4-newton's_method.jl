# newtons_method.jl
# Description:
# Solves various equations using Newton's method.
# Author: adzetto
# Date: 2024-06-18

using Printf

function newtons_method(f, f_prime, x0, tol, max_iter)
    x = x0
    for iter in 1:max_iter
        fx = f(x)
        fpx = f_prime(x)
        if abs(fx) < tol
            return x, iter
        end
        x = x - fx / fpx
    end
    return x, max_iter
end

# Exercise 10: Root near x = -0.5 of e^x - 3x^2 = 0
f1(x) = exp(x) - 3x^2
f1_prime(x) = exp(x) - 6x
root1, iter1 = newtons_method(f1, f1_prime, -0.5, 1e-6, 1000)
@printf("Root near x = -0.5 of e^x - 3x^2 using Newton's method: %.8f (Iterations: %d)\n", root1, iter1)

# Exercise 11: Positive root near x = 4.0 of e^x - 3x^2 = 0
root2, iter2 = newtons_method(f1, f1_prime, 4.0, 1e-6, 1000)
@printf("Positive root near x = 4.0 of e^x - 3x^2 using Newton's method: %.8f (Iterations: %d)\n", root2, iter2)

# Exercise 12: Solve equations from Exercise 4
f2(x) = x^2 - 4
f2_prime(x) = 2x
root3, iter3 = newtons_method(f2, f2_prime, 3, 1e-6, 1000)
@printf("Root of f2 using Newton's method: %.8f (Iterations: %d)\n", root3, iter3)

# Exercise 13a: Derive algorithm for square root of N
function sqrt_newton(N, x0, tol, max_iter)
    x = x0
    for iter in 1:max_iter
        fx = x^2 - N
        fpx = 2x
        if abs(fx) < tol
            return x, iter
        end
        x = (x + N / x) / 2
    end
    return x, max_iter
end

# Example: Find sqrt(16)
root4, iter4 = sqrt_newton(16, 1.0, 1e-6, 1000)
@printf("Square root of 16 using derived Newton's method: %.8f (Iterations: %d)\n", root4, iter4)

# Exercise 13b: Derive similar formulas for third and fourth roots of N
function cbrt_newton(N, x0, tol, max_iter)
    x = x0
    for iter in 1:max_iter
        fx = x^3 - N
        fpx = 3x^2
        if abs(fx) < tol
            return x, iter
        end
        x = (2x + N / x^2) / 3
    end
    return x, max_iter
end

function fourth_root_newton(N, x0, tol, max_iter)
    x = x0
    for iter in 1:max_iter
        fx = x^4 - N
        fpx = 4x^3
        if abs(fx) < tol
            return x, iter
        end
        x = (3x + N / x^3) / 4
    end
    return x, max_iter
end

# Example: Find cbrt(27)
root5, iter5 = cbrt_newton(27, 1.0, 1e-6, 1000)
@printf("Cube root of 27 using derived Newton's method: %.8f (Iterations: %d)\n", root5, iter5)

# Example: Find fourth root of 16
root6, iter6 = fourth_root_newton(16, 1.0, 1e-6, 1000)
@printf("Fourth root of 16 using derived Newton's method: %.8f (Iterations: %d)\n", root6, iter6)

# Exercise 16: Polynomial with a root at x=2 and a triple root at x=1
f3(x) = x^4 - 5x^3 + 9x^2 - 7x + 2
f3_prime(x) = 4x^3 - 15x^2 + 18x - 7

# Starting at x = 2.1
root7a, iter7a = newtons_method(f3, f3_prime, 2.1, 1e-6, 1000)
@printf("Root starting at x = 2.1 using Newton's method: %.8f (Iterations: %d)\n", root7a, iter7a)

# Starting at x = 0.9
root7b, iter7b = newtons_method(f3, f3_prime, 0.9, 1e-6, 1000)
@printf("Root starting at x = 0.9 using Newton's method: %.8f (Iterations: %d)\n", root7b, iter7b)

# Using secant method starting at f(0.9) and f(1.1)
root8, iter8 = secant_method(f3, 0.9, 1.1, 1e-6, 1000)
@printf("Root using secant method starting at x = 0.9 and x = 1.1: %.8f (Iterations: %d)\n", root8, iter8)



#=
Section 1.4
10. Find a root near $x=-0.5$ of $e^x-3 x^2=0$ by Newton's method, to six-digit accuracy.
11. The equation $e^x-3 x^2=0$ has a root not only near $x=-0.5$, but also near $x=4.0$. Find the positive root by Newton's method.
12. Use Newton's method to solve the equations in Exercise 4. How many iterations are required to attain the specified accuracy?
13. a) Use Newton's method on the equation $x^2=N$ to derive the algorithm for the square root of $N$ :
$$
x_{i+1}=\frac{1}{2}\left(x_i+\frac{N}{x_i}\right)
$$
where $x_0$ is an initial approximation to $\sqrt{N}$.
b) Derive similar formulas for the third and fourth roots of $N$.
14. a) If the algorithm of Exercise 13 is applied twice, show that
$$
\sqrt{N}=\frac{A+B}{4}+\frac{N}{A+B}, \quad \text { where } N=A B .
$$
b) Show also that the relative error (error/true value) in (a) is approximately
$$
\frac{1}{8}\left(\frac{A-B}{A+B}\right)^4
$$
15. Expand $f(x)$ about the point $x=a$ in a Taylor series. (See Appendix $\mathrm{A}$ if you have forgotten this.) Using appropriate terms from this, derive the formula for Newton's method.
16. $(x-1)^3(x-2)=x^4-5 x^3+9 x^2-7 x+2=0$ obviously has a root at $x=2$, and a triple root at $x=1$. Beginning with $x=2.1$, use Newton's method once, and observe the degree of improvement. Then start with $x=0.9$, and note the much slower convergence to the triple root even though the initial error is only 0.1 in each case. Use the secant method beginning with $f(0.9)$ and $f(1.1)$, and observe that just one application brings one quite close to the root in contrast to Newton's method. Explain.
=#
