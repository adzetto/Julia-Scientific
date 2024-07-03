module ARMAEstimation

export starma, karma, kalfor, inclu2, regres, forkal

const zero = 0.0
const one = 1.0
const two = 2.0

function starma(ip::Int, iq::Int, ir::Int, np::Int, phi::Vector{Float64}, theta::Vector{Float64}, a::Vector{Float64}, p::Vector{Float64}, v::Vector{Float64}, thetab::Vector{Float64}, xnext::Vector{Float64}, xrow::Vector{Float64}, rbar::Vector{Float64}, nrbar::Int)
    ifault = 0
    if ip < 0
        ifault = 1
    elseif iq < 0
        ifault += 2
    elseif ip == 0 && iq == 0
        ifault = 4
    end
    k = max(iq + 1, ip)
    if ir != k
        ifault = 5
    elseif np != ir * (ir + 1) // 2
        ifault = 6
    elseif nrbar != np * (np - 1) // 2
        ifault = 7
    elseif ir == 1
        ifault = 8
    end
    if ifault != 0
        return ifault
    end

    # Now set A(0), V, and PHI
    for i in 2:ir
        a[i] = zero
        if i > ip
            phi[i] = zero
        end
        v[i] = zero
        if i <= iq + 1
            v[i] = theta[i - 1]
        end
    end
    a[1] = zero
    if ip == 0
        phi[1] = zero
    end
    v[1] = one
    ind = ir
    for j in 2:ir
        vj = v[j]
        for i in j:ir
            ind += 1
            v[ind] = v[i] * vj
        end
    end

    # Now find P(0)
    if ip != 0
        irank = 0
        ssqerr = zero
        rbar .= zero
        p .= zero
        thetab .= zero
        xnext .= zero
        ind = 0
        ind1 = 0
        npr = np - ir
        npr1 = npr + 1
        indj = npr1
        ind2 = npr
        for j in 1:ir
            phij = phi[j]
            xnext[indj] = zero
            indj += 1
            indi = npr1 + j
            for i in j:ir
                ind += 1
                ynext = v[ind]
                phii = phi[i]
                if j != ir
                    xnext[indj] = -phii
                    if i != ir
                        xnext[indi] -= phij
                        ind1 += 1
                        xnext[ind1] = -one
                    end
                end
                xnext[npr1] = -phii * phij
                ind2 += 1
                if ind2 > np
                    ind2 = 1
                end
                xnext[ind2] += one
                inclu2(np, one, xnext, xrow, ynext, p, rbar, thetab, ssqerr, recres, irank)
                xnext[ind2] = zero
                if i != ir
                    xnext[indi] = zero
                    indi += 1
                    xnext[ind1] = zero
                end
            end
        end
        regres(np, nrbar, rbar, thetab, p)

        # Now reorder P
        ind = npr
        for i in 1:ir
            ind += 1
            xnext[i] = p[ind]
        end
        ind = np
        ind1 = npr
        for i in 1:npr
            p[ind] = p[ind1]
            ind -= 1
            ind1 -= 1
        end
        p[1:ir] = xnext[1:ir]
        return ifault
    end

    # P(0) is obtained by backsubstitution for a moving average process
    indn = np + 1
    ind = np + 1
    for i in 1:ir
        for j in 1:i
            ind -= 1
            p[ind] = v[ind]
            if j != 1
                indn -= 1
                p[ind] += p[indn]
            end
        end
    end
    return ifault
end

function karma(ip::Int, iq::Int, ir::Int, phi::Vector{Float64}, theta::Vector{Float64}, a::Vector{Float64}, p::Vector{Float64}, v::Vector{Float64}, n::Int, w::Vector{Float64}, resid::Vector{Float64}, sumlog::Float64, ssq::Float64, iupd::Int, delta::Float64, e::Vector{Float64}, nit::Int)
    ir1 = ir - 1
    e .= zero
    inde = 1

    if nit == 0
        for i in 1:n
            wnext = w[i]

            if iupd != 1 || i != 1
                dt = zero
                if ir != 1
                    dt = p[ir + 1]
                end
                if dt < delta
                    continue
                end
                a1 = a[1]
                if ir != 1
                    for j in 1:ir1
                        a[j] = a[j + 1]
                    end
                end
                a[ir] = zero
                if ip != 0
                    for j in 1:ip
                        a[j] += phi[j] * a1
                    end
                end
                ind = 0
                indn = ir
                for l in 1:ir
                    for j in l:ir
                        ind += 1
                        p[ind] = v[ind]
                        if j != ir
                            indn += 1
                            p[ind] += p[indn]
                        end
                    end
                end
            end

            ft = p[1]
            ut = wnext - a[1]
            if ir != 1
                ind = ir
                for j in 2:ir
                    g = p[j] / ft
                    a[j] += g * ut
                    for l in j:ir
                        ind += 1
                        p[ind] -= g * p[l]
                    end
                end
            end
            a[1] = wnext
            p[1:ir] .= zero
            resid[i] = ut / sqrt(ft)
            e[inde] = resid[i]
            inde += 1
            if inde > iq
                inde = 1
            end
            ssq += ut * ut / ft
            sumlog += log(ft)
        end
        nit = n
        return
    end

    i = 1
    nit = i - 1
    for ii in i:n
        et = w[ii]
        indw = ii
        if ip != 0
            for j in 1:ip
                indw -= 1
                if indw < 1
                    break
                end
                et -= phi[j] * w[indw]
            end
        end
        if iq != 0
            for j in 1:iq
                inde -= 1
                if inde == 0
                    inde = iq
                end
                et -= theta[j] * e[inde]
            end
        end
        e[inde] = et
        resid[ii] = et
        ssq += et * et
        inde += 1
        if inde > iq
            inde = 1
        end
    end
end

function kalfor(m::Int, ip::Int, ir::Int, np::Int, phi::Vector{Float64}, a::Vector{Float64}, p::Vector{Float64}, v::Vector{Float64})
    ir1 = ir - 1
    work = zeros(Float64, ir)
    for l in 1:m
        a1 = a[1]
        if ir != 1
            for i in 1:ir1
                a[i] = a[i + 1]
            end
        end
        a[ir] = zero
        if ip != 0
            for j in 1:ip
                a[j] += phi[j] * a1
            end
        end

        work[1:ir] .= p[1:ir]
        ind = 0
        ind1 = ir
        dt = p[1]
        for j in 1:ir
            phij = phi[j]
            phijdt = phij * dt
            for i in j:ir
                ind += 1
                phii = phi[i]
                p[ind] = v[ind] + phii * phijdt
                if j < ir
                    p[ind] += work[j + 1] * phii
                end
                if i != ir
                    ind1 += 1
                    p[ind] += work[i + 1] * phij + p[ind1]
                end
            end
        end
    end
end

function inclu2(np::Int, weight::Float64, xnext::Vector{Float64}, xrow::Vector{Float64}, ynext::Float64, d::Vector{Float64}, rbar::Vector{Float64}, thetab::Vector{Float64}, ssqerr::Float64, recres::Float64, irank::Int)
    y = ynext
    wt = weight
    xrow[1:np] .= xnext[1:np]
    recres = zero
    ifault = 1
    if wt <= zero
        return
    end
    ifault = 0

    ithisr = 0
    for i in 1:np
        if xrow[i] == zero
            ithisr += np - i
        else
            xi = xrow[i]
            di = d[i]
            dpi = di + wt * xi * xi
            d[i] = dpi
            cbar = di / dpi
            sbar = wt * xi / dpi
            wt = cbar * wt
            if i != np
                i1 = i + 1
                for k in i1:np
                    ithisr += 1
                    xk = xrow[k]
                    rbthis = rbar[ithisr]
                    xrow[k] = xk - xi * rbthis
                    rbar[ithisr] = cbar * rbthis + sbar * xk
                end
            end
            xk = y
            y = xk - xi * thetab[i]
            thetab[i] = cbar * thetab[i] + sbar * xk
            if di == zero
                irank += 1
                return
            end
        end
    end
    ssqerr += wt * y * y
    recres = y * sqrt(wt)
end

function regres(np::Int, nrbar::Int, rbar::Vector{Float64}, thetab::Vector{Float64}, beta::Vector{Float64})
    ithisr = nrbar
    im = np
    for i in 1:np
        bi = thetab[im]
        if im != np
            i1 = i - 1
            jm = np
            for j in 1:i1
                bi -= rbar[ithisr] * beta[jm]
                ithisr -= 1
                jm -= 1
            end
        end
        beta[im] = bi
        im -= 1
    end
end

function forkal(ip::Int, iq::Int, ir::Int, np::Int, ird::Int, irz::Int, id::Int, il::Int, n::Int, nrbar::Int, phi::Vector{Float64}, theta::Vector{Float64}, delta::Vector{Float64}, w::Vector{Float64}, y::Vector{Float64}, amse::Vector{Float64}, a::Vector{Float64}, p::Vector{Float64}, v::Vector{Float64}, resid::Vector{Float64}, e::Vector{Float64}, xnext::Vector{Float64}, xrow::Vector{Float64}, rbar::Vector{Float64}, thetab::Vector{Float64}, store::Vector{Float64})
    ifault = 0
    if ip < 0
        ifault = 1
    elseif iq < 0
        ifault += 2
    elseif ip * ip + iq * iq == 0
        ifault = 4
    end
    k = max(iq + 1, ip)
    if ir != k
        ifault = 5
    elseif np != ir * (ir + 1) // 2
        ifault = 6
    elseif nrbar != np * (np - 1) // 2
        ifault = 7
    elseif id < 0
        ifault = 8
    elseif ird != ir + id
        ifault = 9
    elseif irz != ird * (ird + 1) // 2
        ifault = 10
    elseif il < 1
        ifault = 11
    end
    if ifault != 0
        return ifault
    end

    a[1] = zero
    v[1] = one
    if np == 1
        return ifault
    end
    v[2:np] .= zero
    if iq != 0
        iq1 = iq + 1
        for i in 2:iq1
            v[i] = theta[i - 1]
        end
        for j in 1:iq
            ll = j * (2 * ir + 1 - j) // 2
            for i in j:iq
                lli = ll + i
                v[lli] = theta[i] * theta[j]
            end
        end
    end

    if ir == 1
        p[1] = one / (one - phi[1] * phi[1])
    elseif ir != 1
        starma(ip, iq, ir, np, phi, theta, a, p, v, thetab, xnext, xrow, rbar, nrbar)
    end

    nt = n - id
    if id == 0
        goto170 = true
    else
        for j in 1:id
            nj = n - j
            store[j] = w[nj]
        end
        for i in 1:nt
            aa = zero
            for k in 1:id
                idk = id + i - k
                aa -= delta[k] * w[idk]
            end
            iid = i + id
            w[i] = w[iid] + aa
        end
    end

    if !goto170
        sumlog = zero
        ssq = zero
        iupd = 1
        del = -one
        nit = 0
        karma(ip, iq, ir, phi, theta, a, p, v, nt, w, resid, sumlog, ssq, iupd, del, e, nit)

        sigma = zero
        for j in 1:nt
            sigma += resid[j]^2
        end
        sigma /= nt

        if id != 0
            xrow[1:np] .= p[1:np]
            p[1:irz] .= zero
            ind = 0
            for j in 1:ir
                k = (j - 1) * (id + ir + 1) - (j - 1) * j // 2
                for i in j:ir
                    ind += 1
                    k += 1
                    p[k] = xrow[ind]
                end
            end
            for j in 1:id
                irj = ir + j
                a[irj] = store[j]
            end
        end

        ir2 = ir + 1
        ir1 = ir - 1
        id1 = id - 1
        id2r = 2 * ird
        id2r1 = id2r - 1
        idd1 = 2 * id + 1
        idd2 = idd1 + 1
        i45 = id2r + 1
        idrr1 = ird + 1
        iddr = 2 * id + ir
        jkl = ir * (iddr + 1) // 2
        jkl1 = jkl + 1
        id2r2 = id2r + 2
        ibc = ir * (i45 - ir) // 2

        for l in 1:il
            a1 = a[1]
            if ir != 1
                for i in 1:ir1
                    a[i] = a[i + 1]
                end
            end
            a[ir] = zero
            if ip != 0
                for j in 1:ip
                    a[j] += phi[j] * a1
                end
            end
            if id != 0
                for j in 1:id
                    irj = ir + j
                    a1 += delta[j] * a[irj]
                end
                if id >= 2
                    for i in 1:id1
                        iri1 = ird - i
                        a[iri1 + 1] = a[iri1]
                    end
                end
            end
            a[ir2] = a1

            if id != 0
                for i in 1:id
                    store[i] = zero
                    for j in 1:id
                        ll = max(i, j)
                        k = min(i, j)
                        jj = jkl + (ll - k) + 1 + (k - 1) * (idd2 - k) // 2
                        store[i] += delta[j] * p[jj]
                    end
                end

                for j in 1:id1
                    jj = id - j
                    lk = (jj - 1) * (idd2 - jj) // 2 + jkl
                    lk1 = jj * (idd1 - jj) // 2 + jkl
                    for i in 1:j
                        lk += 1
                        lk1 += 1
                        p[lk1] = p[lk]
                    end
                end
                for j in 1:id1
                    jklj = jkl1 + j
                    irj = ir + j
                    p[jklj] = store[j] + p[irj]
                end
                p[jkl1] = p[1]
                for i in 1:id
                    iri = ir + i
                    p[jkl1] += delta[i] * (store[i] + two * p[iri])
                end
                for i in 1:id
                    iri = ir + i
                    store[i] = p[iri]
                end
                for j in 1:ir
                    kk1 = j * (id2r1 - j) // 2 + ir
                    k1 = (j - 1) * (id2r - j) // 2 + ir
                    for i in 1:id
                        kk = kk1 + i
                        k = k1 + i
                        p[k] = phi[j] * store[i]
                        if j != ir
                            p[k] += p[kk]
                        end
                    end
                end

                for j in 1:ir
                    store[j] = zero
                    kkk = j * (i45 - j) // 2 - id
                    for i in 1:id
                        kkk += 1
                        store[j] += delta[i] * p[kkk]
                    end
                end
                if id >= 2
                    for j in 1:ir
                        k = j * idrr1 - j * (j + 1) // 2 + 1
                        for i in 1:id1
                            k -= 1
                            p[k] = p[k - 1]
                        end
                    end
                end
                for j in 1:ir
                    k = (j - 1) * (id2r - j) // 2 + ir + 1
                    p[k] = store[j] + phi[j] * p[1]
                    if j < ir
                        p[k] += p[j + 1]
                    end
                end
            end
            store[1:ir] = p[1:ir]

            ind = 0
            dt = p[1]
            for j in 1:ir
                phij = phi[j]
                phijdt = phij * dt
                ind2 = (j - 1) * (id2r2 - j) // 2
                ind1 = j * (i45 - j) // 2
                for i in j:ir
                    ind += 1
                    ind2 += 1
                    phii = phi[i]
                    p[ind2] = v[ind] + phii * phijdt
                    if j < ir
                        p[ind2] += store[j + 1] * phii
                    end
                    if i != ir
                        ind1 += 1
                        p[ind2] += store[i + 1] * phij + p[ind1]
                    end
                end
            end

            y[l] = a[1]
            if id != 0
                for j in 1:id
                    irj = ir + j
                    y[l] += a[irj] * delta[j]
                end
            end

            ams = p[1]
            for j in 1:id
                jrj = ibc + (j - 1) * (idd2 - j) // 2
                irj = ir + j
                ams += two * delta[j] * p[irj] + p[jrj + 1] * delta[j]^2
            end
            for j in 1:id1
                j1 = j + 1
                jrk = ibc + 1 + (j - 1) * (idd2 - j) // 2
                for i in j1:id
                    jrj += 1
                    ams += two * delta[i] * delta[j] * p[jrk]
                end
            end
            amse[l] = ams * sigma
        end
    end
end

end # module ARMAEstimation
