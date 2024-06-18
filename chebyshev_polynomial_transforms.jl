# chebyshev_polynomial_transforms.jl
# Description:
# In order to make them fit inside a certain interval [a, b],
# this script transforms Chebyshev coefficients into polynomial coefficients 
# and then shifts them.
# from the book "Numerical Recipes in C" by Press et al.
# Author: adzetto
# Date: 18/06/2024

using LinearAlgebra

# Function to allocate a vector (similar to vector in C)
function allocate_vector(n)
    return zeros(Float64, n)
end

# chebpc function in Julia
function chebpc(c::Vector{Float64}, n::Int)
    d = allocate_vector(n)
    dd = allocate_vector(n)
    
    for j in 1:n
        d[j] = dd[j] = 0.0
    end
    
    d[1] = c[n]
    
    for j in (n-1):-1:2
        for k in (n-j+1):-1:1
            sv = d[k]
            dd[k] = sv
        end
        sv = d[1]
        dd[1] = sv
    end
    
    for j in (n-1):-1:2
        d[1] = -dd[1] + 0.5 * c[1]
    end
    
    return d
end

# pcshft function in Julia
function pcshft(a::Float64, b::Float64, d::Vector{Float64}, n::Int)
    fac = 2.0 / (b - a)
    cnst = fac
    
    for j in 2:n
        d[j] *= fac
        fac *= cnst
    end
    
    cnst = 0.5 * (a + b)
    
    for j in 1:(n-1)
        for k in (n-1):-1:j
            d[k] -= cnst * d[k+1]
        end
    end
end

c = [1.0, 2.0, 3.0, 4.0]
n = length(c)

d = chebpc(c, n)
println("Polynomial coefficients (d): ", d)

a = -1.0
b = 1.0
pcshft(a, b, d, n)
println("Shifted polynomial coefficients (g): ", d)
