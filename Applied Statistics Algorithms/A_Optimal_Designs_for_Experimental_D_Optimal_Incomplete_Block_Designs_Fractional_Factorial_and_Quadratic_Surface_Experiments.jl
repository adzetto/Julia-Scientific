module DOptimalDesign

using LinearAlgebra
using Random

function dopt(x::Matrix{Float64}, dim1::Int, ncand::Int, kin::Int, n::Int, nblock::Int, 
              in::Vector{Int}, blksiz::Vector{Int}, k::Int, rstart::Bool, nrbar::Int,
              d::Vector{Float64}, rbar::Vector{Float64}, picked::Vector{Int}, 
              lndet::Ref{Float64}, xx::Vector{Float64}, tol::Vector{Float64}, 
              zpz::Matrix{Float64}, wk::Vector{Float64}, ifault::Ref{Int})

    one = 1.0
    minus1 = -1.0
    eps = eps(one)
    above1 = 1.0001
    small = 1.0e-4
    hundrd = 100.0

    ifault[] = 0
    if dim1 < ncand ifault[] = 1 end
    if k > n ifault[] += 2 end
    if nrbar < k * (k - 1) // 2 ifault[] += 4 end
    if k != kin + nblock ifault[] += 8 end

    if nblock > 1
        l = 0
        for block in 1:nblock
            l += blksiz[block]
        end
        if n != l ifault[] += 16 end
    else
        if n != blksiz[1] ifault[] += 16 end
    end

    nb = max(1, nblock)
    nin = 0
    for i in 1:nb
        if in[i] < 0 break end
        if in[i] > 0
            if in[i] > blksiz[i] break end
            nin += in[i]
        end
    end
    navail = ncand - nin
    if nin <= n break end
    ifault[] += 32

    if ifault[] != 0 return end

    clear(k, nrbar, d, rbar, ifault)

    tol[1:k] .= eps
    block = 1
    for case in 1:ncand
        getx(x, kin, nblock, block, xx, case)
        for i in 1:k
            tol[i] += abs(xx[i])
        end
    end
    temp = Float64(n) * eps / ncand
    for i in 1:k
        if i <= nblock
            tol[i] = eps
        else
            tol[i] *= temp
        end
    end

    pos = 1
    for block in 1:nb
        if rstart
            last1 = (in[block] + blksiz[block]) // 2
            inc = sqrt(Float64(ncand) + small)
        end
        for i in 1:blksiz[block]
            if rstart && i > in[block]
                rand = rand()
                point = nin + 1 + navail * rand

                if i > last1
                    mxrank = 0
                    lndet[] = -hundrd
                    for cand in nin+1:inc:ncand
                        getx(x, kin, nblock, block, xx, point)
                        modtr2(k, xx, d, rbar, tol, rank, total)
                        if rank < mxrank break end
                        if rank == mxrank && total < lndet[] break end
                        best = point
                        mxrank = rank
                        lndet[] = total * above1
                        point += inc
                        if point > ncand point -= navail end
                    end
                    point = best
                end
                picked[pos] = point
            else
                point = picked[pos]
            end
            getx(x, kin, nblock, block, xx, point)
            modtri(k, one, xx, d, rbar, tol)
            pos += 1
        end
    end

    singm(k, nrbar, d, rbar, tol, wk, ifault)

    if ifault[] == 0 return end

    for pos in 1:k
        if d[pos] < tol[pos] break end
    end

    l = pos - 1
    for i in 1:pos-1
        wk[i] = rbar[l]
        l += k - i - 1
    end
    regcf(k, nrbar, d, rbar, wk, tol, wk, pos-1, ifault)

    bl = 1
    rand = rand()
    case = nin + 1 + navail * rand
    for cand in 1:ncand
        getx(x, kin, nblock, bl, xx, case)
        total = xx[pos] - dot(wk[1:pos-1], xx[1:pos-1])
        if abs(total) > hundrd * tol[pos] break end
        case += 1
        if case > ncand case = nin + 1 end
    end

    if bl == 0
        ifault[] = -1
        return
    end

    bl = 0
    temp = one - small
    pos = in[1] + 1
    for block in 1:nb
        for j in in[block]+1:blksiz[block]
            l = picked[pos]
            getx(x, kin, nblock, block, xx, l)
            bksub2(rbar, k, xx, wk)
            total = dot(wk[1:k], wk[1:k] ./ d[1:k])
            if total < temp
                temp = total
                rempos = pos
                bl = block
            end
            pos += 1
        end
        if block < nblock pos += in[block+1] end
    end

    if bl == 0
        ifault[] = -1
        return
    end

    getx(x, kin, nblock, bl, xx, case)
    modtri(k, one, xx, d, rbar, tol)
    l = picked[rempos]
    getx(x, kin, nblock, bl, xx, l)
    modtri(k, minus1, xx, d, rbar, tol)
    picked[rempos] = case

    for block in 1:nb
        for case in nin+1:ncand
            getx(x, kin, nblock, block, xx, case)
            bksub2(rbar, k, xx, wk)
            temp = dot(wk[1:k], wk[1:k] ./ d[1:k])
            zpz[case, block] = temp
        end
    end

    lastin = 0
    lstout = 0
    change = false
    last = 0
    for block in 1:nb
        first = last + 1 + in[block]
        last += blksiz[block]
        dmax = small
        best = 0

        rand = rand()
        pos = first + (blksiz[block] - in[block]) * rand
        for case in in[block]+1:blksiz[block]
            pos += 1
            if pos > last pos = first end
            i = picked[pos]
            if i == lastin continue end
            getx(x, kin, nblock, block, xx, i)
            bksub2(rbar, k, xx, wk)
            bksub1(rbar, k, wk, wk, tol, d)

            rand = rand()
            j = nin + 1 + navail * rand
            for cand in nin+1:ncand
                j += 1
                if j > ncand j = nin+1 end
                if j == i || j == lstout continue end
                temp = zpz[j, block] - zpz[i, block]
                if temp < dmax continue end
                getx(x, kin, nblock, block, xx, j)
                total = dot(xx[1:k], wk[1:k])
                temp += total^2 - zpz[i, block] * zpz[j, block]
                if temp > dmax
                    dmax = temp * above1
                    best = j
                    rempos = pos
                    drop = i
                end
            end
        end

        if best != 0
            change = true
            if nb == 1
                lastin = best
                lstout = drop
            end

            getx(x, kin, nblock, block, xx, best)
            bksub2(rbar, k, xx, wk)
            bksub1(rbar, k, wk, wk, tol, d)
            modtri(k, one, xx, d, rbar, tol)
            temp = one + zpz[best, block]
            for bl in 1:nb
                for case in 1:ncand
                    getx(x, kin, nblock, bl, xx, case)
                    total = dot(xx[1:k], wk[1:k])
                    zpz[case, bl] -= total^2 / temp
                end
            end

            getx(x, kin, nblock, block, xx, drop)
            bksub2(rbar, k, xx, wk)
            bksub1(rbar, k, wk, wk, tol, d)
            modtri(k, minus1, xx, d, rbar, tol)
            temp = one - zpz[drop, block]
            for bl in 1:nb
                for case in 1:ncand
                    getx(x, kin, nblock, bl, xx, case)
                    total = dot(xx[1:k], wk[1:k])
                    zpz[case, bl] += total^2 / temp
                end
            end
            picked[rempos] = best
        end
    end

    if change return end

    if nblock <= 1 return end

    rpos = nblock * k - nblock * (nblock + 1) // 2 + 1

    last1 = 0
    pos1 = nblock
    dmax = small
    change = false
    for block1 in 1:nblock-1
        first1 = last1 + 1 + in[block1]
        last1 += blksiz[block1]
        last2 = last1
        pos2 = pos1 + k - 1 - block1
        for block2 in block1+1:nblock
            first2 = last2 + 1 + in[block2]
            last2 += blksiz[block2]
            for case1 in in[block1]+1:blksiz[block1]
                posi = first1 - 1 + case1
                i = picked[posi]
                getx(x, kin, nblock, block, xx, i)
                zpz[1:kin, 1] = xx[nblock+1:nblock+kin]
                for case2 in in[block2]+1:blksiz[block2]
                    posj = first2 - 1 + case2
                    j = picked[posj]
                    if i == j continue end
                    getx(x, kin, nblock, block, xx, j)
                    zpz[1:kin, 2] = xx[nblock+1:nblock+kin]
                    delta(kin, zpz[:, 1], zpz[:, 2], rbar[pos1:], rbar[pos2:], 
                          blksiz[block1], blksiz[block2], d[nblock+1:], 
                          rbar[rpos:], temp)
                    if temp > dmax
                        dmax = temp * above1
                        bestb1 = block1
                        bestb2 = block2
                        bestp1 = posi
                        bestp2 = posj
                        change = true
                    end
                end
            end
            pos2 += k - 1 - block2
        end
        pos1 += k - 1 - block1
    end

    if change
        i = picked[bestp1]
        j = picked[bestp2]
        getx(x, kin, nblock, bestb2, xx, i)
        modtri(k, one, xx, d, rbar, tol)
        getx(x, kin, nblock, bestb1, xx, j)
        modtri(k, one, xx, d, rbar, tol)
        getx(x, kin, nblock, bestb1, xx, i)
        modtri(k, minus1, xx, d, rbar, tol)
        getx(x, kin, nblock, bestb2, xx, j)
        modtri(k, minus1, xx, d, rbar, tol)
        picked[bestp1] = j
        picked[bestp2] = i
    end

    lndet[] = sum(log(d[1:k]))
end

function modtri(np::Int, weight::Float64, xrow::Vector{Float64}, d::Vector{Float64}, 
                rbar::Vector{Float64}, tol::Vector{Float64})

    zero = 0.0

    w = weight
    nextr = 1
    for i in 1:np
        if w == zero return end

        xi = xrow[i]
        if abs(xi) < tol[i]
            nextr += np - i
            continue
        end

        di = d[i]
        wxi = w * xi
        dpi = di + wxi * xi

        if dpi < tol[i]
            dpi = zero
            cbar = zero
            sbar = zero
            w = zero
        else
            cbar = di / dpi
            sbar = wxi / dpi
            w = cbar * w
        end

        d[i] = dpi
        for k in i+1:np
            xk = xrow[k]
            xrow[k] -= xi * rbar[nextr]
            rbar[nextr] = cbar * rbar[nextr] + sbar * xk
            nextr += 1
        end
    end
end

function modtr2(np::Int, xrow::Vector{Float64}, d::Vector{Float64}, rbar::Vector{Float64},
                tol::Vector{Float64}, rank::Ref{Int}, lndet::Ref{Float64})

    zero = 0.0
    one = 1.0

    w = one
    rank[] = 0
    lndet[] = zero
    nextr = 1
    for i in 1:np
        if w == zero
            for j in i:np
                if d[j] > tol[j]
                    rank[] += 1
                    lndet[] += log(d[j])
                end
            end
            return
        end

        xi = xrow[i]
        if abs(xi) < tol[i]
            if d[i] > tol[i]
                rank[] += 1
                lndet[] += log(d[i])
            end
            nextr += np - i
            continue
        end

        di = d[i]
        wxi = w * xi
        dpi = di + wxi * xi

        if dpi < tol[i]
            dpi = zero
            cbar = zero
            w = zero
        else
            cbar = di / dpi
            w = cbar * w
            lndet[] += log(dpi)
            rank[] += 1
        end

        for k in i+1:np
            xk = xrow[k]
            xrow[k] -= xi * rbar[nextr]
            nextr += 1
        end
    end
end

function getx(x::Matrix{Float64}, kin::Int, nblock::Int, block::Int, 
              xx::Vector{Float64}, case::Int)
    zero = 0.0
    one = 1.0

    for i in 1:nblock
        if i != block
            xx[i] = zero
        else
            xx[i] = one
        end
    end
    xx[nblock+1:nblock+kin] .= x[case, 1:kin]
end

function bksub1(rbar::Vector{Float64}, k::Int, rhs::Vector{Float64}, 
                soln::Vector{Float64}, tol::Vector{Float64}, d::Vector{Float64})

    zero = 0.0

    pos = k * (k - 1) // 2
    for row in k:-1:1
        if d[row] > tol[row]
            temp = rhs[row] / d[row]
            for col in k:-1:row+1
                temp -= rbar[pos] * soln[col]
                pos -= 1
            end
            soln[row] = temp
        else
            pos -= k - row
            soln[row] = zero
        end
    end
end

function singm(np::Int, nrbar::Int, d::Vector{Float64}, rbar::Vector{Float64}, 
               tol::Vector{Float64}, work::Vector{Float64}, ifault::Ref{Int})

    zero = 0.0

    ifault[] = 0
    if np <= 0 ifault[] = 1 end
    if nrbar < np * (np - 1) // 2 ifault[] += 2 end
    if ifault[] != 0 return end

    work[1:np] .= sqrt.(d[1:np])

    for col in 1:np
        temp = tol[col]
        pos = col - 1
        for row in 1:col-1
            if abs(rbar[pos]) * work[row] < temp rbar[pos] = zero end
            pos += np - row - 1
        end

        if work[col] <= temp
            ifault[] -= 1
            if col < np
                np2 = np - col
                pos2 = pos + np - col + 1
                if np2 > 1
                    modtri(np2, d[col], rbar[pos+1:], d[col+1:], rbar[pos2:], tol)
                else
                    modtri(1, d[col], rbar[pos+1:], d[col+1:], rbar, tol)
                end
                rbar[pos+1:pos2-1] .= zero
            end
            d[col] = zero
        end
    end
end

function xxtr(np::Int, d::Vector{Float64}, rbar::Vector{Float64}, 
              nreq::Int, trace::Ref{Float64}, rinv::Vector{Float64})

    one = 1.0
    zero = 0.0

    inv(np, rbar, nreq, rinv)

    trace[] = zero
    pos = 1
    for row in 1:nreq
        trace[] += one / d[row]
        for col in row+1:nreq
            trace[] += rinv[pos]^2 / d[col]
            pos += 1
        end
    end
end

function inv(np::Int, rbar::Vector{Float64}, nreq::Int, rinv::Vector{Float64})

    zero = 0.0

    pos = nreq * (nreq - 1) // 2
    for row in nreq-1:-1:1
        start = (row - 1) * (np + np - row) // 2 + 1
        for col in nreq:-1:row+1
            pos1 = start
            pos2 = pos
            sum = zero
            for k in row+1:col-1
                pos2 += nreq - k
                sum -= rbar[pos1] * rinv[pos2]
                pos1 += 1
            end
            rinv[pos] = sum - rbar[pos1]
            pos -= 1
        end
    end
end

function clear(np::Int, nrbar::Int, d::Vector{Float64}, rbar::Vector{Float64}, 
               ier::Ref{Int})

    zero = 0.0

    ier[] = 0
    if np < 1 ier[] = 1 end
    if nrbar < np * (np - 1) // 2 ier[] += 2 end
    if ier[] != 0 return end

    d[1:np] .= zero
    rbar[1:nrbar] .= zero
end

function regcf(np::Int, nrbar::Int, d::Vector{Float64}, rbar::Vector{Float64}, 
               thetab::Vector{Float64}, tol::Vector{Float64}, 
               beta::Vector{Float64}, nreq::Int, ier::Ref{Int})

    zero = 0.0

    ier[] = 0
    if np < 1 ier[] = 1 end
    if nrbar < np * (np - 1) // 2 ier[] += 2 end
    if nreq < 1 || nreq > np ier[] += 4 end
    if ier[] != 0 return end

    for i in nreq:-1:1
        if sqrt(d[i]) < tol[i]
            beta[i] = zero
            d[i] = zero
            continue
        end
        beta[i] = thetab[i]
        nextr = (i - 1) * (np + np - i) // 2 + 1
        for j in i+1:nreq
            beta[i] -= rbar[nextr] * beta[j]
            nextr += 1
        end
    end
end

function bksub2(rbar::Vector{Float64}, k::Int, rhs::Vector{Float64}, 
                soln::Vector{Float64})

    for row in 2:k
        temp = rhs[row]
        pos = row - 1
        for col in 1:row-1
            temp -= rbar[pos] * soln[col]
            pos += k - col - 1
        end
        soln[row] = temp
    end
end

function delta(k::Int, xj::Vector{Float64}, xl::Vector{Float64}, 
               xbari::Vector{Float64}, xbark::Vector{Float64}, ni::Int, 
               nk::Int, d::Vector{Float64}, rbar::Vector{Float64}, 
               fn_val::Ref{Float64})

    one = 1.0
    two = 2.0

    const = two - one / ni - one / nk
    for i in 1:k
        temp = xj[i] - xl[i]
        diff[i] = -temp
        a[i] = temp - xbari[i] + xbark[i]
        b[i] = a[i] - const * temp
    end

    bksub2(rbar, k, a, z[:, 1])
    bksub2(rbar, k, b, z[:, 2])
    bksub2(rbar, k, diff, z[:, 3])

    e11 = dot(z[:, 3], z[:, 1] ./ d)
    e12 = dot(z[:, 3], z[:, 3] ./ d)
    e21 = dot(z[:, 2], z[:, 1] ./ d)
    e22 = dot(z[:, 2], z[:, 3] ./ d)

    fn_val[] = (e11 + one) * (e22 + one) - e12 * e21 - one
end

end

using .DOptimalDesign

function ibquad()
    maxf = 6
    maxb = 10
    maxcol = maxf * (maxf + 3) // 2
    kmax = maxcol + maxb
    nrmax = kmax * (kmax - 1) // 2
    mxcand = 3^maxf
    nmax = 100

    blksiz = zeros(Int, maxb)
    in = zeros(Int, maxb)
    picked = zeros(Int, nmax)
    wk = zeros(Float64, kmax)
    d = zeros(Float64, kmax)
    rbar = zeros(Float64, nrmax)
    tol = zeros(Float64, kmax)
    zpz = zeros(Float64, mxcand, maxb)
    design = zeros(Int, nmax)
    lndet = Ref{Float64}(0.0)
    xx = zeros(Float64, kmax)

    println("+++ Incomplete blocks for quad. surfaces +++")
    bel = Char(7)
    
    nfact = 0
    while nfact < 2 || nfact > maxf
        print("Enter no. of factors: ")
        nfact = parse(Int, readline())
        if nfact > maxf || nfact < 2
            println("$bel ** Illegal value entered **")
        end
    end

    print("Enter no. of blocks: ")
    nblock = parse(Int, readline())
    if nblock > maxb || nblock < 0
        println("$bel ** Illegal value entered **")
        return
    end

    print("Enter size of each block: ")
    blksiz = [parse(Int, x) for x in split(readline())]
    n = sum(blksiz)

    print("How many tries?: ")
    nrep = parse(Int, readline())

    Random.seed!(rand(Int))

    kin = nfact * (nfact + 3) // 2
    kfull = kin + nblock
    if n < kfull
        println("$bel ** Design too small to fit model **")
        return
    end
    nrbar = kfull * (kfull - 1) // 2
    ncand = 3^nfact

    open("DESIGN.OPT", "w") do io
        println(io, "Output from DOPT for fitting quadratic surfaces")
        println(io, "$nfact factors at 2 levels")
        if nblock > 1
            println(io, "Block sizes: ", blksiz)
        else
            println(io, "No. of experimental runs = ", n)
        end

        point = 1
        ind = -ones(Int, nfact)
        x = zeros(Float64, mxcand, maxcol)
        
        while true
            x[point, 1:nfact] .= ind[1:nfact]
            pos = nfact + 1
            for i in 1:nfact
                for j in 1:i
                    x[point, pos] = ind[i] * ind[j]
                    pos += 1
                end
            end
            point += 1
            
            i = 1
            while true
                ind[i] += 1
                if ind[i] > 1
                    if i == nfact
                        break
                    end
                    ind[1:i] .= -1
                    i += 1
                else
                    break
                end
            end
            if i == nfact && ind[i] == 2
                break
            end
        end

        dbest = 0.0
        for i in 1:nrep
            ifault = Ref{Int}(0)
            DOptimalDesign.dopt(x, mxcand, ncand, kin, n, nblock, in, blksiz, kfull, 
                                true, nrbar, d, rbar, picked, lndet, xx, tol, zpz, wk, ifault)
            if ifault[] != 0
                println("$bel IFAULT = ", ifault[])
            else
                println(i, " Log Det. = ", lndet[])
            end
            if lndet[] > dbest
                design = copy(picked)
                dbest = lndet[]
            end
        end

        println()
        stdet = exp(dbest - kfull * log(Float64(n)))
        println("Max. log det. = ", dbest, "   Std. det = ", stdet)
        println("Design:")
        println(io, "Max. log det. = ", dbest, "   Std. det = ", stdet)
        println(io, "Design:")
        i2 = 0
        for i in 1:nblock
            i1 = i2 + 1
            i2 += blksiz[i]
            println(io, " BLOCK ", i)
            println(io, "  F A C T O R")
            println(io, "   ", join(1:nfact, " "))
            for j in i1:i2
                println(io, "    ", join(x[design[j], 1:nfact], " "))
            end
        end
        println()
    end
end

function driver()
    maxf = 6
    maxb = 10
    maxcol = maxf * (maxf + 3) // 2
    kmax = maxcol + maxb
    nrmax = kmax * (kmax - 1) // 2
    mxcand = 3^maxf
    nmax = 100

    println("+++ Driver program for DOPT program +++")
    println()

    bel = Char(7)
    
    p = 0
    q = 0
    while p < 0 || q < 0 || (p == 0 && q == 0)
        print("Enter no. of factors with 2 levels: ")
        p = parse(Int, readline())
        if p < 0
            println("$bel ** Number must be positive or zero **")
        end
        print("Enter no. of factors with 3 levels: ")
        q = parse(Int, readline())
        if q < 0
            println("$bel ** Number must be positive or zero **")
        end
        if p == 0 && q == 0
            println("$bel ** Number of factors must be > 0 **")
        end
    end

    nfact = p + q
    print("Enter no. of blocks: ")
    nblock = parse(Int, readline())
    if nblock < 0
        println("$bel ** Negative no. of blocks **")
        return
    end
    nblock = max(1, nblock)

    blksiz = zeros(Int, nblock)
    if nblock == 1
        print("Enter size of the experiment: ")
        blksiz[1] = parse(Int, readline())
    else
        print("Enter size of each block: ")
        blksiz .= [parse(Int, x) for x in split(readline())]
    end
    n = sum(blksiz)

    kin = nfact + q + nfact * (nfact - 1) // 2
    kfull = kin + nblock
    println("No. of parameters in model = $kfull")

    if n < kfull
        println("$bel ** Design too small to fit model **")
        return
    end

    nrbar = kfull * (kfull - 1) // 2
    ncand = (2^p) * (3^q)
    println("No. of candidate points = $ncand")

    print("How many tries?: ")
    nrep = parse(Int, readline())

    Random.seed!(rand(Int))

    open("DESIGN.OPT", "w") do io
        println(io, "Output from DOPT for fractional factorial")
        println(io, "$p factors at 2 levels, and $q factors at 3 levels")
        if nblock > 1
            println(io, "Block sizes: ", blksiz)
        else
            println(io, "No. of experimental runs = ", n)
        end
        println(io, "No. of parameters in full model incl. blocks = $kfull")
        println(io, "No. of candidate points = $ncand")
        println(io, "No. of tries to find optimum design = $nrep")

        point = 1
        ind = -ones(Int, nfact)
        x = zeros(Float64, ncand, kfull)
        
        while true
            x[point, 1:nfact] .= ind[1:nfact]
            pos = nfact
            for i in 1:q
                x[point, pos + i] = ind[i + p]^2
            end
            pos += q
            for i in 1:nfact-1
                for j in i+1:nfact
                    pos += 1
                    x[point, pos] = x[point, i] * x[point, j]
                end
            end
            point += 1
            
            i = 1
            while true
                if i <= p
                    ind[i] = -ind[i]
                    if ind[i] == -1
                        i += 1
                        if i > nfact
                            break
                        end
                    else
                        break
                    end
                else
                    ind[i] += 1
                    if ind[i] > 1
                        ind[i] = -1
                        i += 1
                        if i > nfact
                            break
                        end
                    else
                        break
                    end
                end
            end
            if i == nfact && ind[i] == 2
                break
            end
        end

        d = zeros(Float64, kfull)
        rbar = zeros(Float64, nrbar)
        picked = zeros(Int, n)
        xx = zeros(Float64, kfull)
        tol = zeros(Float64, kfull)
        zpz = zeros(Float64, ncand, nblock)
        wk = zeros(Float64, kfull)
        design = zeros(Int, n)
        lndet = Ref{Float64}(0.0)
        dbest = 0.0
        for i in 1:nrep
            ifault = Ref{Int}(0)
            DOptimalDesign.dopt(x, ncand, ncand, kin, n, nblock, in, blksiz, kfull, 
                                true, nrbar, d, rbar, picked, lndet, xx, tol, zpz, wk, ifault)
            if ifault[] != 0
                println("$bel IFAULT = ", ifault[])
            else
                println(io, "Try no. $i     Log of determinant = ", lndet[])
            end
            if lndet[] > dbest
                design .= copy(picked)
                dbest = lndet[]
            end
        end

        println()
        stdet = exp(dbest - kfull * log(Float64(n)))
        println("Max. log det. = $dbest     Std. det = $stdet")
        println("Design:")
        println(io, "Max. log det. = $dbest     Std. det = $stdet")
        println(io, "Design:")
        i2 = 0
        for i in 1:nblock
            i1 = i2 + 1
            i2 += blksiz[i]
            println(io, " BLOCK $i")
            println(io, "  F A C T O R")
            println(io, "   ", join(1:nfact, " "))
            for j in i1:i2
                println(io, "    ", join(x[design[j], 1:nfact], " "))
            end
        end
        println()
    end
end

ibquad()
driver()
