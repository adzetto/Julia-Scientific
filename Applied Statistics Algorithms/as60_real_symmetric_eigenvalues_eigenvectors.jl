# Algorithm as 60.1 appl.statist. (1973) vol.22 no.2 & Algorithm AS 60.2 appl.statist. (1973) vol.22, no.2
# Translated from FORTRAN by adzetto

module AS60

using LinearAlgebra


const dp = Float64
const zero = 0.0
const one = 1.0
const eta = 1.0e-37
const eps = 1.0e-14

function tdiag(n::Int, a::Matrix{dp})
    """
    Reduces a real symmetric matrix to tridiagonal form.
    
    Arguments:
    n  :: Integer      -- The order of the matrix.
    a  :: Matrix{dp}   -- The symmetric matrix to be reduced.
    
    Returns:
    d  :: Vector{dp}   -- Diagonal elements of the tridiagonal matrix.
    e  :: Vector{dp}   -- Off-diagonal elements of the tridiagonal matrix.
    z  :: Matrix{dp}   -- Orthogonal transformation matrix.
    tol:: Float64      -- Tolerance level used for numerical stability.
    ifault:: Integer   -- Error indicator (0 if no error).
    """
    tol = eta / eps
    ifault = 1

    if n <= 1
        return zeros(dp, n), zeros(dp, n), copy(a), tol, ifault
    end

    ifault = 0
    d = zeros(dp, n)
    e = zeros(dp, n)
    z = copy(a)
    
    for i in n:-1:2
        l = i - 2
        f = z[i, i - 1]
        g = zero
        if l >= 1
            g = sum(z[i, 1:l].^2)
        end
        h = g + f^2

        if g <= tol
            e[i] = f
            d[i] = zero
            continue
        end

        g = sqrt(h)
        if f >= zero
            g = -g
        end
        e[i] = g
        h -= f * g
        z[i, i - 1] = f - g
        f = zero
        for j in 1:l+1
            z[j, i] = z[i, j] / h
            g = zero
            for k in 1:j
                g += z[j, k] * z[i, k]
            end
            if j < l+1
                for k in j+1:l+1
                    g += z[k, j] * z[i, k]
                end
            end
            e[j] = g / h
            f += g * z[j, i]
        end

        hh = f / (h + h)
        for j in 1:l+1
            f = z[i, j]
            g = e[j] - hh * f
            e[j] = g
            for k in 1:j
                z[j, k] -= f * e[k] + g * z[i, k]
            end
        end
        d[i] = h
    end

    d[1] = zero
    e[1] = zero

    for i in 1:n
        l = i - 1
        if d[i] != zero && l != 0
            for j in 1:l
                g = zero
                for k in 1:l
                    g += z[i, k] * z[k, j]
                end
                for k in 1:l
                    z[k, j] -= g * z[k, i]
                end
            end
        end
        d[i] = z[i, i]
        z[i, i] = one
        if l != 0
            for j in 1:l
                z[i, j] = zero
                z[j, i] = zero
            end
        end
    end
    return d, e, z, tol, ifault
end

function lrvt(n::Int, d::Vector{dp}, e::Vector{dp}, z::Matrix{dp})
    """
    Finds eigenvalues and eigenvectors of a tridiagonal matrix using the QL algorithm with implicit shifts.
    
    Arguments:
    n  :: Integer      -- The order of the matrix.
    d  :: Vector{dp}   -- Diagonal elements of the tridiagonal matrix.
    e  :: Vector{dp}   -- Off-diagonal elements of the tridiagonal matrix.
    z  :: Matrix{dp}   -- Orthogonal transformation matrix.
    
    Returns:
    d  :: Vector{dp}   -- Eigenvalues.
    z  :: Matrix{dp}   -- Eigenvectors.
    precis:: Float64   -- Precision level used for numerical stability.
    ifault:: Integer   -- Error indicator (0 if no error).
    """
    mits = 30
    precis = 1.0e-14
    ifault = 2

    if n <= 1
        return d, z, precis, ifault
    end

    ifault = 1
    n1 = n - 1
    e = [e[2:end]; zero]

    b = zero
    f = zero
    for l in 1:n
        jj = 0
        h = precis * (abs(d[l]) + abs(e[l]))
        if b < h
            b = h
        end

        m = l
        for m1 in l:n
            m = m1
            if abs(e[m]) <= b
                break
            end
        end

        if m == l
            d[l] += f
            continue
        end

        while jj < mits
            jj += 1
            p = (d[l + 1] - d[l]) / (2 * e[l])
            r = sqrt(p^2 + one)
            if p < zero
                pr = p - r
            else
                pr = p + r
            end

            h = d[l] - e[l] / pr
            for i in l:n
                d[i] -= h
            end
            f += h

            p = d[m]
            c = one
            s = zero
            for i in m:-1:l+1
                j = i - 1
                g = c * e[j]
                h = c * p
                if abs(p) >= abs(e[j])
                    c = e[j] / p
                    r = sqrt(c^2 + one)
                    e[i] = s * e[j] * r
                    s = one / r
                    c /= r
                else
                    c = p / e[j]
                    r = sqrt(c^2 + one)
                    e[i] = s * p * r
                    s = c / r
                    c = one / r
                end
                p = c * d[j] - s * g
                d[i] = h + s * (c * g + s * d[j])

                for k in 1:n
                    h = z[k, i]
                    z[k, i] = s * z[k, j] + c * h
                    z[k, j] = c * z[k, j] - s * h
                end
            end
            e[l] = s * p
            d[l] = c * p
            if abs(e[l]) <= b
                break
            end
        end

        if jj == mits
            return d, z, precis, ifault
        end

        d[l] += f
    end

    for i in 1:n1
        k = i
        p = d[i]
        for j in i+1:n
            if d[j] > p
                k = j
                p = d[j]
            end
        end

        if k != i
            d[k] = d[i]
            d[i] = p
            for j in 1:n
                p = z[j, i]
                z[j, i] = z[j, k]
                z[j, k] = p
            end
        end
    end

    ifault = 0
    return d, z, precis, ifault
end

end



using .AS60

A = [4.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 4.0]
n = size(A, 1)
d, e, z, tol, ifault_tdiag = AS60.tdiag(n, A)

println("Diagonal elements of the tridiagonal matrix: ", d)
println("Off-diagonal elements of the tridiagonal matrix: ", e)
println("Orthogonal transformation matrix: ", z)
println("Tolerance used: ", tol)
println("Fault code for tdiag: ", ifault_tdiag)

d, z, precis, ifault_lrvt = AS60.lrvt(n, d, e, z)

println("Eigenvalues: ", d)
println("Eigenvectors: ")
println(z)
println("Precision used: ", precis)
println("Fault code for lrvt: ", ifault_lrvt)