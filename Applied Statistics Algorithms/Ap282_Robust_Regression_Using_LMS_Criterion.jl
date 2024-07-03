module RobustRegression

using LinearAlgebra, Random

export mvelms

const dp = Float64

function mvelms(x::Matrix{dp}, y::Vector{dp}, ncas::Int, npre::Int, ioptn::Int, maxtry::Int, nclow::Int, nchigh::Int)
    coeffs = zeros(dp, npre + 1, nchigh - nclow + 1)
    eprmin = fill(Inf, nchigh - nclow + 1)
    resid = zeros(dp, ncas, nchigh - nclow + 1)
    robdsq = zeros(dp, ncas, nchigh - nclow + 1)
    cvemin = fill(Inf, nchigh - nclow + 1)
    ifault = 0

    # Local variables
    itemp = ioptn
    lms = itemp % 2 == 0
    itemp ÷= 2
    mve = itemp % 2 == 0
    itemp ÷= 2
    isint = itemp % 2 == 0
    itemp ÷= 2
    exh = itemp % 2 == 0

    if isint
        nvar = npre + 1
        addcon = 1.0 / nvar
        rpp1 = sqrt(addcon)
        int = 1
    else
        nvar = npre
        addcon = 0.0
        rpp1 = 1.0
        int = 0
    end

    # Check fault parameters
    if nclow < nvar || nclow > nchigh || nchigh > ncas
        ifault = 2
    elseif ncas < nvar
        ifault = 1
    elseif isint && nvar == 1 && mve
        ifault += 4
    end
    if ifault > 0
        return coeffs, eprmin, resid, robdsq, cvemin, ifault
    end

    # Allocate arrays IWORK, DATA & WORK
    iwork = zeros(Int, 4 * nvar + ncas)
    data = zeros(dp, 0:nvar, nvar + ncas)
    work = zeros(dp, 3 * ncas + nvar)

    exact = false
    j2 = ncas
    j3 = 2 * ncas
    j4 = 3 * ncas
    i2 = nvar
    i3 = 2 * nvar
    i4 = 3 * nvar
    i5 = i4 + ncas
    nvr2 = 2 * nvar
    deter = rpp1
    powmed = npre * 0.5
    for i in 1:nvar
        iwork[i] = ncas + i - nvar
    end
    iptr = 1
    ncrang = nchigh - nclow + 1
    eprmin .= Inf
    cvemin .= Inf
    lfresh = 0
    ncount = 0

    # Target variable, if present, is treated in standardized form so that exact fit can be detected more easily
    if lms
        ixlo = 0
        for i in 1:ncas
            work[i] = y[i]
        end

        # Sort entries
        sort!(work[1:ncas])
        ymed = (work[ncas ÷ 2 + 1] + work[(ncas + 1) ÷ 2]) * 0.5
        for i in 1:ncas
            work[i] = abs(y[i] - ymed)
        end
        sort!(work[1:ncas])
        ymad = (work[ncas ÷ 2 + 1] + work[(ncas + 1) ÷ 2]) * 0.5
        if ymad == 0.0
            ymad = 1.0
        end
    else
        ixlo = 1
    end
    detadj = 0.0
    for j in 1:npre
        for i in 1:ncas
            work[i] = abs(x[i, j])
        end
        sort!(work[1:ncas])
        work[j4 + j + int] = (work[ncas ÷ 2 + 1] + work[(ncas + 1) ÷ 2]) * 0.5
        if work[j4 + j + int] == 0.0
            work[j4 + j + int] = work[ncas]
        end
        detadj += log10(work[j4 + j + int])
    end
    if isint
        work[j4 + 1] = 1.0
    end

    # Data is transferred to work area; initial simplex tableau is set up
    refresh!(data, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, i2, i3, i5, 0, isint, int, lms, deter, lfresh)

    # Initial basis is set up. Checks are made that initial basis members are in general position
    if exh
        iwork[i2 + 1] = 1
        for i in 1:nvar
            j = iwork[i2 + i]
            for ii in i:nvar
                if abs(data[ii, j]) >= 1e-9
                    pivot!(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
                    iwork[i3 + i] = j
                    if ii != i
                        swap!(data, nvar, ncas, ii, i)
                    end
                    break
                end
            end
            if j == ncas
                ifault += 8
                if ifault >= 24
                    ifault -= 16
                end
                return coeffs, eprmin, resid, robdsq, cvemin, ifault
            else
                ifault = 16
                iwork[i2 + i] += 1
                continue
            end
        end
        iptr = nvar
    else
        for i in 1:ncas
            iwork[i4 + i] = i
        end
        shuffle!(iwork[i4 + 1:i4 + ncas])
        inxptr = 0
        for i in 1:nvar
            inxptr += 1
            j = iwork[i4 + inxptr]
            for ii in i:nvar
                if abs(data[ii, j]) >= 1e-9
                    pivot!(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
                    iwork[i3 + i] = j
                    if ii != i
                        swap!(data, nvar, ncas, ii, i)
                    end
                    iwork[i2 + i] = j
                    break
                end
            end
            if inxptr == ncas
                ifault += 8
                if ifault >= 24
                    ifault -= 16
                end
                return coeffs, eprmin, resid, robdsq, cvemin, ifault
            else
                ifault = 16
                continue
            end
        end
    end

    while true
        if exh
            iwork[i2 + iptr] += 1
            if iwork[i2 + iptr] > iwork[iptr]
                iptr -= 1
                if iptr == 0
                    break
                else
                    continue
                end
            else
                if abs(data[iptr, iwork[i2 + iptr]]) >= 1e-9
                    ncount += 1
                    lfresh += 1
                    pivot!(data, ixlo, nvar, ncas, iwork, i3, iptr, iwork[i2 + iptr], deter)
                    if lfresh > 50
                        deter = rpp1
                        refresh!(data, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, i2, i3, i5, nvar, isint, int, lms, deter, lfresh)
                    end
                    if exh && iptr < nvar
                        iptr += 1
                        iwork[i2 + iptr] = iwork[i2 + iptr - 1]
                        continue
                    end
                else
                    if lfresh > nvr2
                        deter = rpp1
                        refresh!(data, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, i2, i3, i5, iptr - 1, isint, int, lms, deter, lfresh)
                        if abs(data[iptr, iwork[i2 + iptr]]) >= 1e-9
                            continue
                        end
                        ifault = 16
                        if iptr == nvar
                            continue
                        else
                            deter = rpp1
                            refresh!(data, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, i2, i3, i5, iptr - 1, isint, int, lms, deter, lfresh)
                            j = iwork[i2 + iptr]
                            for ii in iptr:nvar
                                if abs(data[ii, j]) >= 1e-9
                                    pivot!(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
                                    iwork[i3 + iptr] = j
                                    if ii != iptr
                                        swap!(data, nvar, ncas, ii, iptr)
                                    end
                                    break
                                end
                            end
                            iwork[i2 + iptr] += 1
                            if iwork[i2 + iptr] > iwork[iptr]
                                iptr -= 1
                                if iptr == 0
                                    break
                                else
                                    iwork[i2 + iptr] += 1
                                end
                            end
                            continue
                        end
                    else
                        ifault = 16
                        if iptr == nvar
                            continue
                        else
                            deter = rpp1
                            refresh!(data, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, i2, i3, i5, iptr - 1, isint, int, lms, deter, lfresh)
                            j = iwork[i2 + iptr]
                            for ii in iptr:nvar
                                if abs(data[ii, j]) >= 1e-9
                                    pivot!(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
                                    iwork[i3 + iptr] = j
                                    if ii != iptr
                                        swap!(data, nvar, ncas, ii, iptr)
                                    end
                                    break
                                end
                            end
                            iwork[i2 + iptr] += 1
                            if iwork[i2 + iptr] > iwork[iptr]
                                iptr -= 1
                                if iptr == 0
                                    break
                                else
                                    iwork[i2 + iptr] += 1
                                end
                            end
                            continue
                        end
                    end
                end
            end
        else
            if ncount > maxtry
                break
            end
            iptr = mod(iptr, nvar) + 1
            inxptr += 1
            if inxptr > ncas
                for i in 1:nvar
                    iwork[i4 + ncas - nvar + i] = iwork[i4 + i]
                    iwork[i4 + i] = iwork[i2 + i]
                end
                shuffle!(iwork[i4 + 1:i4 + ncas])
                inxptr = nvar + 1
            end
            ncol = iwork[i4 + inxptr]
            if abs(data[iptr, ncol]) < 1e-9
                if lfresh > nvr2
                    deter = rpp1
                    refresh!(data, ixlo, x, ncas, npre, y, ymad, nvar, work, j4, iwork, i2, i3, i5, iptr - 1, isint, int, lms, deter, lfresh)
                    if abs(data[iptr, ncol]) < 1e-9
                        ifault = 16
                        continue
                    end
                else
                    ifault = 16
                    continue
                end
            end
            iwork[i2 + iptr] = ncol
        end
    end

    if mve
        # Convert robust distances and scaling back to original scale
        for k in 1:ncrang
            cvemin[k] += detadj
            for j in 1:ncas
                if robdsq[j, k] > -Inf
                    robdsq[j, k] = 10.0 ^ robdsq[j, k]
                else
                    robdsq[j, k] = 0.0
                end
            end
        end
    end

    if lms
        for k in 1:ncrang
            eprmin[k] *= ymad
            for i in 1:npre
                coeffs[i + int, k] /= work[i + j4 + int]
            end
        end
    end

    return coeffs, eprmin, resid, robdsq, cvemin, ifault
end

function refresh!(data::Matrix{dp}, ixlo::Int, x::Matrix{dp}, ncas::Int, npre::Int, y::Vector{dp}, ymad::dp, nvar::Int, work::Vector{dp}, j4::Int, iwork::Vector{Int}, i2::Int, i3::Int, i5::Int, iup::Int, isint::Bool, int::Int, lms::Bool, deter::dp, lfresh::Int)
    # Subroutine to refresh first iup entries of simplex basis
    lfresh = 0
    if isint
        data[1, 1:ncas] .= 1.0
    end

    for i in 1:ncas
        data[int + 1:int + npre, i] .= x[i, 1:npre] ./ work[j4 + int + 1:j4 + int + npre]
    end

    for i in 1:nvar
        data[1:nvar, ncas + i] .= 0.0
        data[i, ncas + i] = 1.0
    end

    if lms
        data[0, 1:ncas] .= y[1:ncas] ./ ymad
        data[0, ncas + 1:ncas + nvar] .= 0.0
    end

    for i in 1:iup
        j = iwork[i2 + i]
        for ii in i:nvar
            if abs(data[ii, j]) >= 1e-9
                pivot!(data, ixlo, nvar, ncas, iwork, i5, ii, j, deter)
                if ii != i
                    swap!(data, nvar, ncas, ii, i)
                end
                iwork[i3 + i] = j
                break
            end
        end
    end

    return
end

function pivot!(x::Matrix{dp}, ixlo::Int, nord::Int, ncas::Int, iwork::Vector{Int}, i3::Int, nrow::Int, ncol::Int, deter::dp)
    # Subroutine to pivot entry corresponding to (nrow, ncol) into simplex tableau
    pivt = x[nrow, ncol]
    deter *= pivt
    nhigh = nord + ncas
    for j in 1:nhigh
        if j != ncol
            fmult = x[nrow, j] / pivt
            for i in ixlo:nord
                if i != nrow
                    x[i, j] -= fmult * x[i, ncol]
                end
            end
            x[nrow, j] = fmult
        end
    end
    x[ixlo:nord, ncol] .= 0.0
    x[nrow, ncol] = 1.0
    iwork[i3 + nrow] = ncol

    return
end

function perm!(index::Vector{Int}, iaa::Int, n::Int, iab::Int)
    # Subroutine to return pseudorandom permutation of n - iab entries in vector index starting at iaa + 1
    m = n - iab

    while m > 1
        j = rand(1:m)
        # Swap entries in index corresponding to j and m, offset by iab
        if j != m
            index[iaa + j + iab], index[iaa + m + iab] = index[iaa + m + iab], index[iaa + j + iab]
        end
        m -= 1
    end

    return
end

function swap!(data::Matrix{dp}, nvar::Int, ncas::Int, ir1::Int, ir2::Int)
    # Subroutine to swap rows ir1 and ir2 of matrix data
    nhigh = nvar + ncas
    for i in 1:nhigh
        data[ir1, i], data[ir2, i] = data[ir2, i], data[ir1, i]
    end

    return
end

function sortsub!(ra::Vector{dp}, n::Int, ilow::Int)
    # Subroutine to sort entries ilow + 1 through ilow + n in ascending order from a vector ra using heapsort
    l = n ÷ 2 + 1
    ir = n

    while true
        if l > 1
            l -= 1
            rra = ra[l + ilow]
        else
            rra = ra[ir + ilow]
            ra[ir + ilow] = ra[1 + ilow]
            ir -= 1
            if ir == 1
                ra[1 + ilow] = rra
                break
            end
        end

        i = l
        j = l + l

        while j <= ir
            if j < ir
                if ra[j + ilow] < ra[j + 1 + ilow]
                    j += 1
                end
            end
            if rra < ra[j + ilow]
                ra[i + ilow] = ra[j + ilow]
                i = j
                j += j
            else
                j = ir + 1
            end
        end

        ra[i + ilow] = rra
    end

    return
end

end  # module RobustRegression

using .RobustRegression

function main()
    ncas = 50
    npre = 5
    ioptn = 1
    maxtry = 100
    nclow = 5
    nchigh = 20

    x = rand(RobustRegression.dp, ncas, npre)
    y = rand(RobustRegression.dp, ncas)

    coeffs, eprmin, resid, robdsq, cvemin, ifault = RobustRegression.mvelms(x, y, ncas, npre, ioptn, maxtry, nclow, nchigh)

    println("Coefficients:")
    println(coeffs)

    println("Eprmin:")
    println(eprmin)

    println("Residuals:")
    println(resid)

    println("Robdsq:")
    println(robdsq)

    println("Cvemin:")
    println(cvemin)

    println("Ifault:")
    println(ifault)
end

main()
