# Filename: error_calculations.jl
# Description: Julia code to calculate absolute and relative errors for given functions.
# Author: adzetto
# Date: 2024-06-28

# Function to calculate absolute and relative errors for a given function
function calculate_errors(f, df, x, Δx)
    # Calculate the function value at x
    fx = f(x)
    # Calculate the derivative value at x
    dfx = df(x)
    
    # Calculate the absolute error
    Δf = abs(dfx) * Δx
    # Calculate the relative error
    relative_error = Δf / abs(fx)
    
    return Δf, relative_error
end

# Part (i): Relative error for v = (u1^p * u2^q * u3^r) / (u4^s * u5^t)
function relative_error_v(u, Δu, p, q, r, s, t)
    relative_error = abs(p) * (Δu[1] / u[1]) + abs(q) * (Δu[2] / u[2]) +
                     abs(r) * (Δu[3] / u[3]) + abs(s) * (Δu[4] / u[4]) +
                     abs(t) * (Δu[5] / u[5])
    return relative_error
end

# Part (ii): Relative error for f(x) = 2x^5 - 3x + 2 at x = 1 with Δx = 0.005
f1(x) = 2x^5 - 3x + 2
df1(x) = 10x^4 - 3
x1 = 1.0
Δx1 = 0.005
Δf1, relative_error_f1 = calculate_errors(f1, df1, x1, Δx1)
println("Part (ii): Absolute error in f(x) = 2x^5 - 3x + 2 at x = $x1: $Δf1")
println("Part (ii): Relative error in f(x) = 2x^5 - 3x + 2 at x = $x1: $relative_error_f1")

# Part (iii): Absolute and relative errors for y = (1.42x + 3.45) / (x + 0.75) at x = 0.5 ± 0.1
f2(x) = (1.42x + 3.45) / (x + 0.75)
df2(x) = (1.42 * (x + 0.75) - (1.42x + 3.45)) / (x + 0.75)^2
x2 = 0.5
Δx2 = 0.1
Δf2, relative_error_f2 = calculate_errors(f2, df2, x2, Δx2)
println("Part (iii): Absolute error in y = (1.42x + 3.45) / (x + 0.75) at x = $x2: $Δf2")
println("Part (iii): Relative error in y = (1.42x + 3.45) / (x + 0.75) at x = $x2: $relative_error_f2")

# Example for Part (i)
u = [1.0, 2.0, 3.0, 4.0, 5.0]
Δu = [0.01, 0.02, 0.03, 0.04, 0.05]
p, q, r, s, t = 2, 3, 4, 5, 6
relative_error = relative_error_v(u, Δu, p, q, r, s, t)
println("Part (i): Relative error for v = (u1^p * u2^q * u3^r) / (u4^s * u5^t): $relative_error")
