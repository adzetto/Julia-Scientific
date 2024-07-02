# Algorithm AS7, Applied Statistics, vol.17, 1968, p.198.
# N.B. Argument W has been removed.
# Forms in c( ) as lower triangle, a generalised inverse
# of the positive semi-definite symmetric matrix a( )
# order n, stored as lower triangle
# Translated from FORTRAN by adzetto


module SymmetricMatrixInversion

using LinearAlgebra


#     chol(a::Vector{Float64}, n::Int, nn::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int})
# 
# This function computes the Cholesky decomposition of a symmetric positive definite matrix `a`. The Cholesky decomposition is a factorization of `a` into a lower triangular matrix `u` and its conjugate transpose `u'`. The function modifies the input matrix `a` in-place to store the lower triangular matrix `u`. The function also returns the number of null elements in `u` in the `nullty` argument and an error code in the `ifault` argument.
# 
# # Arguments
# - `a::Vector{Float64}`: The input matrix `a` stored as a one-dimensional array in column-major order.
# - `n::Int`: The order of the matrix `a`.
# - `nn::Int`: The length of the array `a`.
# - `u::Vector{Float64}`: The output lower triangular matrix `u` stored as a one-dimensional array in column-major order.
# - `nullty::Ref{Int}`: A reference to an integer that will store the number of null elements in `u`.
# - `ifault::Ref{Int}`: A reference to an integer that will store the error code.
# 
# # Returns
# - `nothing`: If an error occurs, the function returns `nothing`.
# - `nullty::Int`: The number of null elements in `u`.
# - `ifault::Int`: An error code. If `ifault` is 0, the function executed successfully.
# 

function chol(a::Vector{Float64}, n::Int, nn::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int})
    eta = 1e-9
    eta2 = eta * eta
    zero = 0.0

    ifault[] = 1
    if n <= 0
        return
    end

    ifault[] = 3
    if nn < div(n * (n + 1), 2)
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


#     syminv(a::Vector{Float64}, n::Int, nn::Int, c::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int})
# 
# This function computes the inverse of a symmetric positive definite matrix `a`. The inverse is computed using the Cholesky decomposition method and is stored in-place in the input matrix `c`. The function also returns the number of null elements in `c` in the `nullty` argument and an error code in the `ifault` argument.
# 
# # Arguments
# - `a::Vector{Float64}`: The input matrix `a` stored as a one-dimensional array in column-major order.
# - `n::Int`: The order of the matrix `a`.
# - `nn::Int`: The length of the array `a`.
# - `c::Vector{Float64}`: The output matrix `c` (the inverse of `a`) stored as a one-dimensional array in column-major order.
# - `nullty::Ref{Int}`: A reference to an integer that will store the number of null elements in `c`.
# - `ifault::Ref{Int}`: A reference to an integer that will store the error code.
# 
# # Returns
# - `nothing`: If an error occurs, the function returns `nothing`.
# - `nullty::Int`: The number of null elements in `c`.
# - `ifault::Int`: An error code. If `ifault` is 0, the function executed successfully.
# 
# # Description
# The function first computes the Cholesky decomposition of `a` and then uses this decomposition to compute the inverse of `a`, which is stored in `c`. The computation is done in-place, and the original content of `a` is overwritten by its Cholesky decomposition. The function handles errors through the `ifault` reference argument.

function syminv(a::Vector{Float64}, n::Int, nn::Int, c::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int})
    w = Vector{Float64}(undef, n)
    zero = 0.0
    one = 1.0

    chol(a, n, nn, c, nullty, ifault)
    if ifault[] != 0
        return
    end

    ndiag = nn
    for irow in n:-1:1
        l = ndiag
        if c[ndiag] == zero
            for j in irow:n
                c[l] = zero
                l += j
            end
            continue
        end

        for i in irow:n
            w[i] = c[l]
            l += i
        end

        icol = n
        jcol = nn
        mdiag = nn
        while true
            l = jcol
            x = zero
            if icol == irow
                x = one / w[irow]
            end
            k = n
            while true
                if k == irow
                    break
                end
                x -= w[k] * c[l]
                k -= 1
                l -= 1
                if l < 1
                    l = mdiag - k
                end
            end

            c[l] = x / w[irow]
            if icol == irow
                break
            end

            mdiag -= icol
            icol -= 1
            jcol -= 1
        end

        ndiag -= irow
    end

    ifault[] = 0
end


#     subchl(a::Vector{Float64}, b::Vector{Int}, n::Int, u::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int}, det::Ref{Float64}=Ref(1.0))
# 
# This function performs a modified Cholesky decomposition on a symmetric matrix `a`, represented in a compact form. The decomposition is used to compute the determinant of the matrix and to check for its positive definiteness. The function modifies the input matrix `a` in-place to store the Cholesky factors. It also returns the determinant of `a` in the `det` argument, the number of null elements in the Cholesky factor in the `nullty` argument, and an error code in the `ifault` argument.
# 
# # Arguments
# - `a::Vector{Float64}`: The input matrix `a` stored as a one-dimensional array.
# - `b::Vector{Int}`: An array that helps in mapping the compact storage of `a` to its full matrix form.
# - `n::Int`: The order of the matrix `a`.
# - `u::Vector{Float64}`: The output Cholesky factor stored as a one-dimensional array.
# - `nullty::Ref{Int}`: A reference to an integer that will store the number of null elements in the Cholesky factor.
# - `ifault::Ref{Int}`: A reference to an integer that will store the error code.
# - `det::Ref{Float64}`: A reference to a float that will store the determinant of the matrix `a`.
# 
# # Returns
# - `nothing`: If an error occurs or the matrix is not positive definite, the function returns `nothing`.
# 
# # Description
# The function first checks for the positive definiteness of the matrix `a` and computes its Cholesky decomposition. The decomposition is stored in `u`, and the determinant is calculated as the product of the squares of the diagonal elements of `u`. The function sets `ifault` to 0 if the decomposition is successful, to 1 if `n` is non-positive, and to 2 if the matrix is not positive definite. The number of null elements in the Cholesky factor is returned in `nullty`.

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
        ij = div(b[icol] * (b[icol] - 1), 2)
        ii = ij + b[icol]
        x = eta2 * a[ii]
        l = 0
        w = 0.0
        kk = 0.0

        for irow in 1:icol
            kk = div(b[irow] * (b[irow] + 1), 2)
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


#     subinv(a::Vector{Float64}, nm::Int, b::Vector{Int}, n::Int, c::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int}, det::Ref{Float64}=Ref(1.0))
# 
# This function computes the inverse of a symmetric matrix `a` using a modified Cholesky decomposition. The matrix `a` is represented in a compact form, and its inverse is stored in the array `c`. The function also calculates the determinant of `a` and checks for its positive definiteness.
# 
# # Arguments
# - `a::Vector{Float64}`: The input matrix `a` stored as a one-dimensional array.
# - `nm::Int`: The maximum order of the matrix that can be processed.
# - `b::Vector{Int}`: An array that helps in mapping the compact storage of `a` to its full matrix form.
# - `n::Int`: The order of the matrix `a`.
# - `c::Vector{Float64}`: The output array where the inverse of `a` will be stored.
# - `nullty::Ref{Int}`: A reference to an integer that will store the number of null elements in the Cholesky factor.
# - `ifault::Ref{Int}`: A reference to an integer that will store the error code.
# - `det::Ref{Float64}`: A reference to a float that will store the determinant of the matrix `a`.
# 
# # Returns
# - `nothing`: If an error occurs, the function returns `nothing`.
# 
# # Description
# The function first checks the input parameters for validity. It then calls `subchl` to perform a modified Cholesky decomposition of `a`. If the decomposition is successful, the function proceeds to compute the inverse of `a` by back substitution and stores it in `c`. The determinant of `a` is also calculated if requested. The function sets `ifault` to different values based on the type of error encountered.

function subinv(a::Vector{Float64}, nm::Int, b::Vector{Int}, n::Int, c::Vector{Float64}, nullty::Ref{Int}, ifault::Ref{Int}, det::Ref{Float64}=Ref(1.0))
    w = Vector{Float64}(undef, n)
    zero = 0.0
    one = 1.0

    ifault[] = 3
    if n > nm || b[1] < 1 || b[1] > nm - n + 1
        return
    end

    if n > 1
        for i in 2:n
            if b[i] <= b[i-1] || b[i] > nm - n + i
                return
            end
        end
    end

    nrow = n
    ifault[] = 1
    if nrow <= 0
        return
    end

    ifault[] = 0

    if isassigned(det)
        subchl(a, b, nrow, c, nullty, ifault, det)
    else
        subchl(a, b, nrow, c, nullty, ifault)
    end
    if ifault[] != 0
        return
    end

    nn = nrow * (nrow + 1) // 2
    ndiag = nn
    for irow in nrow:-1:1
        if c[ndiag] == zero
            l = ndiag
            for j in irow:nrow
                c[l] = zero
                l += j
            end
            continue
        end

        l = ndiag
        for i in irow:nrow
            w[i] = c[l]
            l += i
        end

        icol = nrow
        jcol = nn
        mdiag = nn
        while true
            l = jcol
            x = zero
            if icol == irow
                x = one / w[irow]
            end
            k = nrow
            while true
                if k == irow
                    break
                end
                x -= w[k] * c[l]
                k -= 1
                l -= 1
                if l < 1
                    l = mdiag - k
                end
            end

            c[l] = x / w[irow]
            if icol == irow
                break
            end
            mdiag -= icol
            icol -= 1
            jcol -= 1
        end
        ndiag -= irow
    end
end

end


using .SymmetricMatrixInversion

function test_as6_and_as7()
    n = 3
    nn = div(n * (n + 1), 2)
    a = rand(Float64, nn)
    c = Vector{Float64}(undef, nn)
    aa = Vector{Float64}(undef, nn)
    u = Vector{Float64}(undef, nn)
    b = collect(1:n)
    nullty = Ref{Int}(0)
    ifault = Ref{Int}(0)
    det = Ref{Float64}(1.0)

    # Generate random lower-triangular matrix.
    c = rand(Float64, nn)

    # aa = c * c'
    det[] = 1.0
    end_row = nn
    destn = nn
    for row in n:-1:1
        pos1 = end_row
        for col in row:-1:1
            pos2 = destn
            if col == row
                det[] *= c[destn]
            end
            total = c[pos1] * c[pos2]
            for k in col-1:-1:1
                pos1 -= 1
                pos2 -= 1
                total += c[pos1] * c[pos2]
            end
            aa[destn] = total
            destn -= 1
            pos1 -= 1
        end
        end_row = destn
    end

    println("Determinant of original Cholesky factor = ", det[])
    println("Expected determinant of inv(A'A)        = ", 1.0 / det[]^2)

    a = aa

    # Form the Cholesky factorization, then invert & calculate inv(A).
    SymmetricMatrixInversion.syminv(a, n, nn, u, nullty, ifault)
    if nullty[] != 0
        println("NULLTY = ", nullty[])
    end
    if ifault[] != 0
        println("IFAULT = ", ifault[])
    end

    # Reverse the process
    b = collect(1:n)

    # Apart from signs, u should be identical to the original c matrix.
    SymmetricMatrixInversion.subinv(u, n, b, n, a, nullty, ifault, det)
    if nullty[] != 0
        println("NULLTY = ", nullty[])
    end
    if ifault[] != 0
        println("IFAULT = ", ifault[])
    end
    println("Determinant now = ", det[])

    # Calculate maximum absolute difference between elements of A & AA.
    maxerr = maximum(abs.(a .- aa))
    println("Max. abs. error = ", maxerr)
end

test_as6_and_as7()
