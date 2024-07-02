# Date: 2024-07-02  Time: 10:12:00
# ALGORITHM AS 6  APPL. STATIST. VOL.17, NO.1
# Translated from FORTRAN by Muhammet Yağcıoğlu

module CholeskyFactorization

# chol(a::Vector{Float64}, n::Int, nn::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int})
#
# Given a symmetric matrix order `n` as lower triangle in `a`, calculates an upper triangle `u` such that `u' * u = a`.
# The input matrix `a` must be positive semi-definite. The factor `eta` determines the effective zero for pivot.
#
# Arguments
# - `a`: Input, positive definite matrix stored in lower-triangular form.
# - `n`: Input, the order of matrix `a`.
# - `nn`: Input, the size of `a` and `u` arrays, must be >= n*(n+1)/2.
# - `u`: Output, a lower triangular matrix such that `u * u' = a`. `a` and `u` may occupy the same locations.
# - `nullty`: Output, the rank deficiency of `a`.
# - `ifault`: Output, error indicator.
#               - 1 if n < 1.
#               - 2 if `a` is not positive semi-definite.
#               - 3 if `nn` < n*(n+1)/2.
#               - 0 otherwise.

function chol(a::Vector{Float64}, n::Int, nn::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int})
    eta = 1e-9
    eta2 = eta * eta
    zero = 0.0

    ifault[] = 1
    if n <= 0
        return
    end

    ifault[] = 3
    if nn < n * (n + 1) // 2
        return
    end

    ifault[] = 2
    nullty[] = 0
    j = 1
    k = 0
    ii = 0

    for icol in 1:n
        ii += icol
        x = eta2 * a[ii]
        l = 0
        kk = 0
        w = 0.0
        for irow in 1:icol
            kk += irow
            k += 1
            w = a[k]
            m = j

            for i in 1:irow
                l += 1
                if i == irow
                    break
                end
                w -= u[l] * u[m]
                m += 1
            end

            if irow == icol
                break
            end

            if u[l] != zero
                u[k] = w / u[l]
            else
                if w * w > abs(x * a[kk])
                    return
                end
                u[k] = zero
            end
        end

        if abs(w) > abs(eta * a[k])
            if w < zero
                return
            end
            u[k] = sqrt(w)
        else
            u[k] = zero
            nullty[] += 1
        end
        j += icol
    end

    ifault[] = 0
end

# Given a symmetric matrix of order `n` as lower triangle in `a`, calculates an upper triangle `u` such that `u' * u` is the sub-matrix
# of `a` whose rows and columns are specified in the integer array `b`. `u` may coincide with `a`. The input matrix `a` must be positive semi-definite.
#
# Arguments
# - `a`: Input, positive definite matrix stored in lower-triangular form.
# - `b`: Input, integer array specifying rows and columns.
# - `n`: Input, the order of matrix `a`.
# - `u`: Output, a lower triangular matrix such that `u * u' = a`. `a` and `u` may occupy the same locations.
# - `nullty`: Output, the rank deficiency of `a`.
# - `ifault`: Output, error indicator.
#               - 1 if n <= 0.
#               - 2 if `a` is not positive semi-definite.
#               - 0 otherwise.
# - `det`: Optional output, determinant of the matrix.

function subchl(a::Vector{Float64}, b::Vector{Int}, n::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int}, det::Ref{Float64}=Ref(1.0))
    eta = 1e-14
    eta2 = eta * eta
    zero = 0.0
    one = 1.0
    ifault[] = 1
    if n <= 0
        return
    end
    ifault[] = 2
    nullty[] = 0
    det[] = one
    j = 1
    k = 0
    for icol in 1:n
        ij = b[icol] * (b[icol] - 1) // 2
        ii = ij + b[icol]
        x = eta2 * a[ii]
        l = 0
        for irow in 1:icol
            kk = b[irow] * (b[irow] + 1) // 2
            k += 1
            jj = ij + b[irow]
            w = a[jj]
            m = j
            for i in 1:irow
                l += 1
                if i == irow
                    break
                end
                w -= u[l] * u[m]
                m += 1
            end
            if irow == icol
                break
            end
            if u[l] != zero
                u[k] = w / u[l]
            else
                if w * w > abs(x * a[kk])
                    return
                end
                u[k] = zero
            end
        end
        if abs(w) > abs(eta * a[kk])
            if w < zero
                return
            end
            u[k] = sqrt(w)
        else
            u[k] = zero
            nullty[] += 1
        end
        j += icol
        det[] *= u[k] * u[k]
    end
    ifault[] = 0
end


function subchl(a::Vector{Float64}, b::Vector{Int}, n::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int}, det::Ref{Float64}=Ref(1.0))
    eta = 1e-14
    eta2 = eta * eta
    zero = 0.0
    one = 1.0

    ifault[] = 1
    if n <= 0
        return
    end

    ifault[] = 2
    nullty[] = 0
    det[] = one

    j = 1
    k = 0

    for icol in 1:n
        ij = b[icol] * (b[icol] - 1) // 2
        ii = ij + b[icol]
        x = eta2 * a[ii]
        l = 0

        for irow in 1:icol
            kk = b[irow] * (b[irow] + 1) // 2
            k += 1
            jj = ij + b[irow]
            w = a[jj]
            m = j

            for i in 1:irow
                l += 1
                if i == irow
                    break
                end
                w -= u[l] * u[m]
                m += 1
            end

            if irow == icol
                break
            end

            if u[l] != zero
                u[k] = w / u[l]
            else
                if w * w > abs(x * a[kk])
                    return
                end
                u[k] = zero
            end
        end

        if abs(w) > abs(eta * a[kk])
            if w < zero
                return
            end
            u[k] = sqrt(w)
        else
            u[k] = zero
            nullty[] += 1
        end

        j += icol
        det[] *= u[k] * u[k]
    end

    ifault[] = 0
end

end

using .CholeskyFactorization

a = [4.0, 1.0, 1.0, 2.0, 1.0, 2.0]
n = 3
nn = length(a)
u = zeros(Float64, nn)
nullty = Ref(0)
ifault = Ref(0)

CholeskyFactorization.chol(a, n, nn, u, nullty, ifault)
println("Cholesky factor u: ", u)
println("Nullity: ", nullty[])
println("Ifault: ", ifault[])