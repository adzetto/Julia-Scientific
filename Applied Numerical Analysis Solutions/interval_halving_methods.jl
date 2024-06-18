# interval_halving_methods.jl
# Description:
# Solves various equations using the interval halving (bisection) method.
# Author: adzetto
# Date: 2024-06-18

using Printf

function bisection_method(f, a, b, tol, max_iter)
    for iter in 1:max_iter
        c = (a + b) / 2
        if f(c) == 0 || (b - a) / 2 < tol
            return c, iter
        end
        if sign(f(c)) == sign(f(a))
            a = c
        else
            b = c
        end
    end
    return (a + b) / 2, max_iter
end

# 1. Finding the root of e^x - 3x using interval halving.
f1(x) = exp(x) - 3x
root1, iter1 = bisection_method(f1, 0, 1, 1e-6, 6)
@printf "Root of e^x - 3x after 6 iterations: %.8f (Iterations: %d)\n" root1, iter1

# Find the number of iterations required for 4 and 8 significant figures.
tol_4sig = 0.00005
tol_8sig = 0.000000005
_, iter_4sig = bisection_method(f1, 0, 1, tol_4sig, 1000)
_, iter_8sig = bisection_method(f1, 0, 1, tol_8sig, 1000)
@printf "Iterations needed for 4 significant places: %d\n" iter_4sig
@printf "Iterations needed for 8 significant places: %d\n" iter_8sig

# 2. Finding zeros of (x-0.4)(x-0.6).
f2(x) = x^2 - x + 0.24
root2a, _ = bisection_method(f2, 0.3, 0.5, 1e-6, 1000)
root2b, _ = bisection_method(f2, 0.5, 0.7, 1e-6, 1000)
@printf "Root near 0.4: %.8f\n" root2a
@printf "Root near 0.6: %.8f\n" root2b

# Error bound and actual error after 5 iterations from interval [0.5, 1.0].
root2c, iter5 = bisection_method(f2, 0.5, 1.0, 1e-6, 5)
error_bound_5 = (1.0 - 0.5) / (2^5)
actual_error_5 = abs(root2c - 0.6)
@printf "Bound to error after 5 iterations: %.8f\n" error_bound_5
@printf "Actual error after 5 iterations: %.8f\n" actual_error_5

# 3. Intersection of y = x - 2 and y = ln(x).
f3(x) = log(x) - x + 2
root3, _ = bisection_method(f3, 1, 3, 1e-4, 1000)
@printf "Intersection of y = x - 2 and y = ln(x): %.8f\n" root3

# 4a. Smallest positive root of e^x - x - 2.
f4a(x) = exp(x) - x - 2
root4a, _ = bisection_method(f4a, 0, 2, 1e-6, 1000)
@printf "Smallest positive root of e^x - x - 2: %.8f\n" root4a

# 4b. Smallest positive root of x^3 - x^2 - 2x + 1.
f4b(x) = x^3 - x^2 - 2x + 1
root4b, _ = bisection_method(f4b, 0, 1, 1e-6, 1000)
@printf "Smallest positive root of x^3 - x^2 - 2x + 1: %.8f\n" root4b

# 4c. Smallest positive root of 2e^(-x) - sin(x).
f4c(x) = 2exp(-x) - sin(x)
root4c, _ = bisection_method(f4c, 0, 2, 1e-6, 1000)
@printf "Smallest positive root of 2e^(-x) - sin(x): %.8f\n" root4c

# 4d. Smallest positive root of 3x^3 + 4x^2 - 8x - 1.
f4d(x) = 3x^3 + 4x^2 - 8x - 1
root4d, _ = bisection_method(f4d, 0, 1, 1e-6, 1000)
@printf "Smallest positive root of 3x^3 + 4x^2 - 8x - 1: %.8f\n" root4d



#=
Section 1.2
1. The equation $e^x-3 x$ has a root at $r=0.61906129$. Beginning with the interval $[0,1]$, use six iterations of the method of halving the interval to find this root. How many iterations would you need to evaluate the root correct to four significant places-that is, $|x-r|<$ 0.00005 ? How many for eight places?
2. The quadratic $(x-0.4)(x-0.6)=x^2-x+0.24$ has zeros at $x=0.4$ and $x=0.6$, of course. Observe that the endpoints of the interval $[0,1]$ are not satisfactory to begin the interval-halving method. Graph the function, and from this deduce the boundaries of intervals that will converge to each of the zeros. If the endpoints of the interval $[0.5,1.0]$ are used to begin the search, what is a bound to the error after five iterations? What is the actual error after five repetitions of interval halving?
3. Interval halving applies to any continuous function, not just to polynomials. Find where the graphs of $y=x-2$ and $y=\ln x$ intersect by finding the root of $\ln x-x+2=0$ correct to four decimals.
4. Use interval halving to find the smallest positive root of these equations. In each case first determine a suitable interval, then compute the root with relative accuracy of $0.5 \%$.
a) $e^x-x-2=0$
b) $x^3-x^2-2 x+1=0$
c) $2 e^{-x}-\sin x=0$
d) $3 x^3+4 x^2-8 x-1=0$
=#