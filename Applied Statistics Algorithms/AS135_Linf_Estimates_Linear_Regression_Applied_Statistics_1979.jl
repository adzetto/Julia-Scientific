"""
    lfnorm(n::Int, m::Int, x::Matrix{Float64}, y::Vector{Float64})

Compute Min-Max (L-infinity) estimates for linear multiple regression.

This function implements Algorithm AS 135 from Applied Statistics (1979) Vol.28, No. 1.

# Arguments
- `n::Int`: Number of cases (observations)
- `m::Int`: Number of predictors (independent variables)
- `x::Matrix{Float64}`: Values of the predictor variables, one case per row
- `y::Vector{Float64}`: Values of the dependent variable

# Returns
- `beta::Vector{Float64}`: Final estimates of the regression coefficients
- `z::Float64`: The least maximum absolute deviation
- `ky::Int`: The number of iterations
- `ifault::Int`: Error indicator (0 for normal termination, 1 if the observation matrix has less than full rank)

# Note
This implementation is based on the FORTRAN 90 version by Alan Miller.
Latest revision of the original code: 5 January 1999
"""
function lfnorm(n::Int, m::Int, x::Matrix{Float64}, y::Vector{Float64})
    # Initialize output variables
    beta = zeros(Float64, m)
    z = 0.0
    ky = 0
    ifault = 0

    # Initialize working variables
    hilo = zeros(Float64, m)
    xrxf = zeros(Float64, m)
    xsxf = zeros(Float64, m)
    lu = zeros(Float64, m, m)
    ibase = zeros(Int, m)
    indx = collect(1:m)

    # Set numerical precision constants
    acu = sqrt(eps(Float64))
    big = sqrt(prevfloat(Inf64))

    # Main algorithm implementation
    m1 = m - 1
    intl = true
    kkk = 1

    # Set up initial LU decomposition
    update!(kkk, x, lu, ibase, indx, intl, n, m, ifault)
    ifault != 0 && return beta, z, ky, ifault
    intl = false
    irow = kkk

    # Calculate initial beta values
    for ii in 1:m
        k = indx[ii]
        k1 = ibase[ii]
        beta[k] = ii == 1 ? y[k1] / lu[k, 1] : (y[k1] - sum(lu[indx[1:ii-1], ii] .* beta[indx[1:ii-1]])) / lu[k, ii]
    end
    for ii in 1:m1
        k1 = m - ii
        k = indx[k1]
        beta[k] -= sum(lu[indx[m-ii+1:m], k1] .* beta[indx[m-ii+1:m]])
    end

    # Main iterative loop
    while true
        # Search for and set first violated constraint
        irow += 1
        irow > n && return beta, z, ky, ifault
        dev1 = sum(x[irow, :] .* beta) - y[irow]
        abs(dev1) < acu && continue
        sigr = sign(dev1)
        rrr = irow

        # Adjust for the R-th constraint
        for ii in 1:m
            k = indx[ii]
            xrxf[ii] = x[rrr, k] - (ii > 1 ? sum(lu[indx[1:ii-1], ii] .* xrxf[1:ii-1]) : 0.0)
            ii == m && (xrxf[ii] /= lu[k, ii])
        end
        hilo[m] = -sign(sigr * xrxf[m])
        sumxr = sigr - hilo[m] * xrxf[m]
        for ii in 1:m1
            k1 = m - ii
            k = indx[k1]
            xrxf[k1] = (xrxf[k1] - sum(lu[k, m-ii+1:m] .* xrxf[m-ii+1:m])) / lu[k, k1]
            hilo[k1] = -sign(sigr * xrxf[k1])
            sumxr -= hilo[k1] * xrxf[k1]
        end
        z = abs(dev1 / sumxr)

        # Start of main iterative loop
        while true
            sss = 0
            deviat = acu

            # Calculate beta values
            for ii in 1:m
                k = indx[ii]
                k1 = ibase[ii]
                beta[k] = ii == 1 ? (y[k1] + z * hilo[1]) / lu[k, 1] :
                          (y[k1] + z * hilo[ii] - sum(lu[indx[1:ii-1], ii] .* beta[indx[1:ii-1]])) / lu[k, ii]
            end
            for ii in 1:m1
                k1 = m - ii
                k = indx[k1]
                beta[k] -= sum(lu[indx[m-ii+1:m], k1] .* beta[indx[m-ii+1:m]])
            end

            # Calculate residuals and find most violated S-th constraint
            for i in 1:n
                yest = sum(x[i, :] .* beta)
                dev1 = abs(y[i] - yest) - z
                if dev1 > deviat
                    ydev = yest - y[i]
                    deviat = dev1
                    sss = i
                end
            end

        # Check if at optimum
        sss == 0 && return beta, z, ky, ifault

        # Set up information on the S-th constraint
        sigs = sign(sum(x[sss, :] .* beta) - y[sss])
        for ii in 1:m
            k = indx[ii]
            xsxf[ii] = x[sss, k] - (ii > 1 ? sum(lu[indx[1:ii-1], ii] .* xsxf[1:ii-1]) : 0.0)
            ii == m && (xsxf[ii] /= lu[k, m])
        end
        sumxs = -sigs + hilo[m] * xsxf[m]
        for ii in 1:m1
            k1 = m - ii
            k = indx[k1]
            xsxf[k1] = (xsxf[k1] - sum(lu[k, m-ii+1:m] .* xsxf[m-ii+1:m])) / lu[k, k1]
            sumxs += hilo[k1] * xsxf[k1]
        end

            # Search for minimum ratio
            while true
                kkk = 0
                ratio = big
                for i in 1:m
                    if sigs * sign(xsxf[i]) != hilo[i] || abs(xsxf[i]) < acu
                        continue
                    end
                    test = abs(xrxf[i] / xsxf[i])
                    if test < ratio
                        ratio = test
                        kkk = i
                    end
                end

                # Process the movement of constraints
                if kkk == 0
                    delta = abs(deviat / sumxs)
                    div = abs(sumxr) - 2.0
                    if div < acu
                        sigr, sigs = sigs, -sigr
                        xrxf, xsxf = xsxf, xrxf
                        sumxr, sumxs = -sumxs, -sumxr + sigr + sigr
                        rrr, sss = sss, rrr
                        deviat = abs(sumxs * delta) - 2.0 * z
                        z += delta
                    else
                        swing = 2.0 * z / div
                        if swing < delta
                            sigr, sigs = sigs, -sigr
                            xrxf, xsxf = xsxf, xrxf
                            sumxr, sumxs = -sumxs, -sumxr + sigr + sigr
                            rrr, sss = sss, rrr
                            deviat = abs(sumxs * delta) - 2.0 * z
                            z += delta
                        else
                            sigr = sigs
                            xrxf = copy(xsxf)
                            sumxr = -sumxs
                            z += delta
                            rrr = sss
                            break
                        end
                    end
                else
                    delta = abs(xrxf[kkk] * deviat / (xrxf[kkk] * sumxs + xsxf[kkk] * sumxr))
                    top = -2.0 * z * xrxf[kkk]
                    div = xrxf[kkk]^2 + hilo[kkk] * sumxr
                    if sign(top) != sign(div) || abs(div) < acu
                        hilo[kkk] = sigs
                        sumxr = sigr
                        xrxf[kkk] /= xsxf[kkk]
                        sumxr -= hilo[kkk] * xrxf[kkk]
                        for i in 1:m
                            if i != kkk
                                xrxf[i] -= xsxf[i] * xrxf[kkk]
                                sumxr -= hilo[i] * xrxf[i]
                            end
                        end
                        ibase[kkk] = sss
                        update!(kkk, x, lu, ibase, indx, intl, n, m, ifault)
                        ifault != 0 && return beta, z, ky, ifault
                        z += delta
                        ky += 1
                        break
                    else
                        swing = top / div
                        if swing < delta
                            z += swing
                            deviat -= swing * abs(sumxs + xsxf[kkk] * sumxr / xrxf[kkk])
                            sumxr += 2.0 * hilo[kkk] * xrxf[kkk]
                            sumxs -= 2.0 * hilo[kkk] * xsxf[kkk]
                            hilo[kkk] = -hilo[kkk]
                        else
                            hilo[kkk] = sigs
                            sumxr = sigr
                            xrxf[kkk] /= xsxf[kkk]
                            sumxr -= hilo[kkk] * xrxf[kkk]
                            for i in 1:m
                                if i != kkk
                                    xrxf[i] -= xsxf[i] * xrxf[kkk]
                                    sumxr -= hilo[i] * xrxf[i]
                                end
                            end
                            ibase[kkk] = sss
                            update!(kkk, x, lu, ibase, indx, intl, n, m, ifault)
                            ifault != 0 && return beta, z, ky, ifault
                            z += delta
                            ky += 1
                            break
                        end
                    end
                end
            end
        end
    end
end

"""
    update!(kkk::Int, x::Matrix{Float64}, lu::Matrix{Float64}, ibase::Vector{Int}, 
            indx::Vector{Int}, intl::Bool, n::Int, m::Int, ifault::Int)

Update LU decomposition matrix.

This function implements Algorithm AS 135.1 from Applied Statistics (1979) Vol.28, No. 1.

# Arguments
- `kkk::Int`: Current iteration index
- `x::Matrix{Float64}`: Input matrix of predictor variables
- `lu::Matrix{Float64}`: LU decomposition matrix (modified in-place)
- `ibase::Vector{Int}`: Base indices (modified in-place)
- `indx::Vector{Int}`: Index vector (modified in-place)
- `intl::Bool`: Flag for initial iteration
- `n::Int`: Number of cases (observations)
- `m::Int`: Number of predictors (independent variables)
- `ifault::Int`: Error indicator (modified in-place)

# Note
This function modifies its arguments in-place.
"""
function update!(kkk::Int, x::Matrix{Float64}, lu::Matrix{Float64}, ibase::Vector{Int}, 
                 indx::Vector{Int}, intl::Bool, n::Int, m::Int, ifault::Int)
    acu = sqrt(eps(Float64))
    irow = 0

    for ii in kkk:m
        if !intl
            irow = ibase[ii]
        else
            irow += 1
            if irow > n
                ifault = 1
                return
            end
        end

        lu[:, ii] = x[irow, :]

        if ii > 1
            for icol in 1:ii-1
                k = indx[icol]
                subt = lu[k, ii]
                for i in icol+1:m
                    k = indx[i]
                    lu[k, ii] -= subt * lu[k, icol]
                end
            end
        end

        pivot = acu
        kk = 0
        for i in ii:m
            k = indx[i]
            if abs(lu[k, ii]) > pivot
                pivot = abs(lu[k, ii])
                kk = i
            end
        end

        if kk == 0
            if intl
                continue
            else
                ifault = 1
                return
            end
        end

        indx[kk], indx[ii] = indx[ii], indx[kk]

        if intl
            ibase[ii] = irow
        end

        if ii < m
            for i in ii+1:m
                k = indx[i]
                lu[k, ii] /= lu[indx[ii], ii]
            end
        end
    end

    kkk = irow
end

# Test program for lfnorm function
function test_lfnorm()
    println("Enter number of predictor variables: ")
    m = parse(Int, readline())
    
    if m < 1
        println("** Must be greater than zero **")
        return
    end

    n = 3 * m

    # Generate artificial data with coefficients 1, 2, ..., m
    x = [ones(n) rand(n, m-1)]
    true_beta = collect(1.0:m)
    y = x * true_beta + rand(n) .- 0.5

    beta, z, niter, ier = lfnorm(n, m, x, y)

    if ier > 0
        println("** X-matrix has rank < m **")
    else
        println("Regression coefficients:")
        println(beta)
        println("No. of iterations = ", niter)
        println("Least max. abs. deviation = ", z)
    end
end

# Run the test program
test_lfnorm()