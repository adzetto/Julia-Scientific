module MinMaxRegression

export lfnorm

using LinearAlgebra

function lfnorm(n::Int, m::Int, x::Array{Float64, 2}, y::Vector{Float64})
    # Outputs
    beta = zeros(Float64, m)
    z = 0.0
    ky = 0
    ifault = 0

    # Constants
    acu = sqrt(eps(Float64))
    big = sqrt(typemax(Float64))
    zero = 0.0
    one = 1.0
    two = 2.0

    # Local variables
    hilo = zeros(Float64, m)
    xrxf = zeros(Float64, m)
    xsxf = zeros(Float64, m)
    lu = zeros(Float64, m, m)
    indx = collect(1:m)
    ibase = zeros(Int, m)

    m1 = m - 1
    intl = true
    kkk = 1

    # Subroutine update call
    ifault = update!(kkk, x, lu, ibase, indx, intl, n, m)
    if ifault != 0
        return beta, z, ky, ifault
    end

    intl = false
    irow = kkk

    # Calculate beta value
    k = indx[1]
    k1 = ibase[1]
    beta[k] = y[k1] / lu[k, 1]
    for ii in 2:m
        k = indx[ii]
        k1 = ibase[ii]
        beta[k] = y[k1]
        ii1 = ii - 1
        for i in 1:ii1
            kk = indx[i]
            beta[k] -= lu[kk, ii] * beta[kk]
        end
        beta[k] /= lu[k, ii]
    end
    for ii in 1:m1
        k1 = m - ii
        k = indx[k1]
        for i in 1:ii
            kk = m - i + 1
            k2 = indx[kk]
            beta[k] -= lu[k2, k1] * beta[k2]
        end
    end

    # Search for and set first violated constraint as R-th constraint
    while true
        irow += 1
        if irow > n
            return beta, z, ky, ifault
        end

        dev1 = zero
        for i in 1:m
            dev1 += x[irow, i] * beta[i]
        end
        dev1 -= y[irow]
        if abs(dev1) < acu
            continue
        end

        sigr = sign(one, dev1)
        rrr = irow

        # Adjust for the R-th constraint
        k = indx[1]
        xrxf[1] = x[rrr, k]
        for ii in 2:m
            k = indx[ii]
            xrxf[ii] = x[rrr, k]
            ii1 = ii - 1
            for i in 1:ii1
                xrxf[ii] -= lu[indx[i], ii] * xrxf[i]
            end
        end
        k = indx[m]
        xrxf[m] /= lu[k, m]
        hilo[m] = sign(one, -sigr * xrxf[m])
        sumxr = sigr - hilo[m] * xrxf[m]
        for ii in 1:m1
            k1 = m - ii
            k = indx[k1]
            for i in 1:ii
                k2 = m - i + 1
                xrxf[k1] -= lu[k, k2] * xrxf[k2]
            end
            xrxf[k1] /= lu[k, k1]
            hilo[k1] = sign(one, -sigr * xrxf[k1])
            sumxr -= hilo[k1] * xrxf[k1]
        end
        z = abs(dev1 / sumxr)

        # Main iterative loop begins
        while true
            sss = 0
            deviat = acu

            # Calculate beta value
            k = indx[1]
            k1 = ibase[1]
            beta[k] = (y[k1] + z * hilo[1]) / lu[k, 1]
            for ii in 2:m
                k = indx[ii]
                k1 = ibase[ii]
                beta[k] = y[k1] + z * hilo[ii]
                ii1 = ii - 1
                for i in 1:ii1
                    kk = indx[i]
                    beta[k] -= lu[kk, ii] * beta[kk]
                end
                beta[k] /= lu[k, ii]
            end
            for ii in 1:m1
                k1 = m - ii
                k = indx[k1]
                for i in 1:ii
                    kk = m - i + 1
                    k2 = indx[kk]
                    beta[k] -= lu[k2, k1] * beta[k2]
                end
            end

            # Calculate residuals
            for i in 1:n
                yest = dot(x[i, 1:m], beta[1:m])
                dev1 = abs(y[i] - yest) - z
                if dev1 <= deviat
                    continue
                end
                ydev = yest - y[i]
                deviat = dev1
                sss = i
            end

            # Check if at optimum
            if sss == 0
                return beta, z, ky, ifault
            end

            # Set up information on the S-th constraint
            sigs = sign(one, ydev)
            k = indx[1]
            xsxf[1] = x[sss, k]
            for ii in 2:m
                k = indx[ii]
                xsxf[ii] = x[sss, k]
                ii1 = ii - 1
                for i in 1:ii1
                    xsxf[ii] -= lu[indx[i], ii] * xsxf[i]
                end
            end
            k = indx[m]
            xsxf[m] /= lu[k, m]
            sumxs = -sigs + hilo[m] * xsxf[m]
            for ii in 1:m1
                k1 = m - ii
                k = indx[k1]
                for i in 1:ii
                    k2 = m - i + 1
                    xsxf[k1] -= lu[k, k2] * xsxf[k2]
                end
                xsxf[k1] /= lu[k, k1]
                sumxs += hilo[k1] * xsxf[k1]
            end

            # Search for minimum ratio
            kkk = 0
            ratio = big
            for i in 1:m
                if sigs * sign(one, xsxf[i]) != hilo[i] || abs(xsxf[i]) < acu
                    continue
                end
                test = abs(xrxf[i] / xsxf[i])
                if test >= ratio
                    continue
                end
                ratio = test
                kkk = i
            end

            # Check if R-th constraint moves interior
            if kkk == 0
                delta = abs(deviat / sumxs)

                # Calculate the largest tolerable delta
                div = abs(sumxr) - two
                if div < acu
                    continue
                end
                swing = two * z / div
                if swing >= delta
                    continue
                end

                # Switch R and S constraint indicators
                save = sumxs
                sumxs = -sumxr + sigr + sigr
                sumxr = -save
                save = sigr
                sigr = sigs
                sigs = -save
                deviat = abs(sumxs * delta) - two * z
                z += delta
                xsxf, xrxf = xrxf, xsxf
                rrr, sss = sss, rrr
                continue
            end

            # Process the movement of the K-th constraint
            delta = abs(xrxf[kkk] * deviat / (xrxf[kkk] * sumxs + xsxf[kkk] * sumxr))
            top = -two * z * xrxf[kkk]
            div = xrxf[kkk] * xrxf[kkk] + hilo[kkk] * sumxr
            if sign(one, top) != sign(one, div) || abs(div) < acu
                continue
            end
            swing = top / div

            # Check to see if the K-th constraint swings across
            if swing >= delta
                continue
            end
            z += swing
            deviat -= swing * abs(sumxs + xsxf[kkk] * sumxr / xrxf[kkk])
            sumxr += two * hilo[kkk] * xrxf[kkk]
            sumxs -= two * hilo[kkk] * xsxf[kkk]
            hilo[kkk] = -hilo[kkk]
            continue

            # Update XRXF and the LU of the current basis
            hilo[kkk] = sigs
            sumxr = sigr
            xrxf[kkk] /= xsxf[kkk]
            sumxr -= hilo[kkk] * xrxf[kkk]
            for i in 1:m
                if i == kkk
                    continue
                end
                xrxf[i] -= xsxf[i] * xrxf[kkk]
                sumxr -= hilo[i] * xrxf[i]
            end
            ibase[kkk] = sss

            # Update LU decomposition
            ifault = update!(kkk, x, lu, ibase, indx, intl, n, m)
            if ifault != 0
                return beta, z, ky, ifault
            end
            z += delta
            ky += 1
        end
    end

    return beta, z, ky, ifault
end

function update!(kkk::Int, x::Array{Float64, 2}, lu::Array{Float64, 2}, ibase::Vector{Int}, indx::Vector{Int}, intl::Bool, n::Int, m::Int)
    # Local variables
    acu = sqrt(eps(Float64))
    ifault = 0

    irow = 0
    for ii in kkk:m
        if !intl
            irow = ibase[ii]
        else
            irow += 1
            if irow > n
                ifault = 1
                return ifault
            end
        end
        lu[1:m, ii] = x[irow, 1:m]

        # Set up representation of incoming row
        if ii != 1
            ii1 = ii - 1
            for icol in 1:ii1
                k = indx[icol]
                subt = lu[k, ii]
                for j in icol + 1:m
                    k = indx[j]
                    lu[k, ii] -= subt * lu[k, icol]
                end
            end
        end

        # Find maximum entry
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
            continue
        end

        # Switch order
        isave = indx[kk]
        indx[kk] = indx[ii]
        indx[ii] = isave

        # Put into columns of LU one at a time
        if intl
            ibase[ii] = irow
        end
        if ii != m
            for i in ii + 1:m
                k = indx[i]
                lu[k, ii] /= lu[isave, ii]
            end
        end
    end
    kkk = irow

    return ifault
end

end # module MinMaxRegression

using .MinMaxRegression

# Example usage
n = 9
m = 3
x = [1.0 1.0 1.0; 1.0 2.0 4.0; 1.0 3.0 9.0; 1.0 4.0 16.0; 1.0 5.0 25.0; 1.0 6.0 36.0; 1.0 7.0 49.0; 1.0 8.0 64.0; 1.0 9.0 81.0]
y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

beta, z, ky, ifault = MinMaxRegression.lfnorm(n, m, x, y)

println("Regression coefficients: ", beta)
println("Least max. abs. deviation: ", z)
println("Number of iterations: ", ky)
println("Ifault: ", ifault)
