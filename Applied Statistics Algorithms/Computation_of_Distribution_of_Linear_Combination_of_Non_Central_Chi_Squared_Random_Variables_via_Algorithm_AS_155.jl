module NoncChisq
using Printf

const dp = Float64
const zero = 0.0
const one = 1.0
const two = 2.0
const pi = 3.14159265358979
const aln28 = 0.0866

mutable struct QFCommon
    aintl::Float64
    ersm::Float64
    sigsq::Float64
    almax::Float64
    almin::Float64
    amean::Float64
    c::Float64
    icount::Int
    ir::Int
    lim::Int
    ndtsrt::Bool
    fail::Bool
    QFCommon() = new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, true, false)
end

const qfcom = QFCommon()

function iround(x::Float64)
    half = 0.5
    return Int(x + sign(x) * half)
end

function qf(alb::Vector{Float64}, anc::Vector{Float64}, n::Vector{Int}, irr::Int, sigma::Float64, cc::Float64, lim1::Int, acc::Float64, ith::Vector{Int}, trace::Vector{Float64}, ifault::Ref{Int})
    local fn_val::Float64
    local j, nj, nt, ntm
    local acc1, almx, xlim, xnt, xntm, utx, tausq, sd, aintv, aintv1, x, up, un, d1, d2, alj, ancj
    
    qfcom.c = cc
    qfcom.ir = irr
    qfcom.lim = lim1
    trace .= zero
    ifault[] = 0
    qfcom.icount = 0
    qfcom.aintl = zero
    qfcom.ersm = zero
    fn_val = -one
    acc1 = acc
    qfcom.ndtsrt = true
    qfcom.fail = false

    xlim = lim1
    qfcom.sigsq = sigma^2
    sd = qfcom.sigsq
    qfcom.almax = zero
    qfcom.almin = zero
    qfcom.amean = zero
    j = 1

    while j <= irr
        nj = n[j]
        alj = alb[j]
        ancj = anc[j]
        if nj < 0 || ancj < zero
            ifault[] = 3
            return fn_val
        end

        sd += alj^2 * (2 * nj + 4 * ancj)
        qfcom.amean += alj * (nj + ancj)
        if qfcom.almax < alj
            qfcom.almax = alj
        end
        if qfcom.almin > alj
            qfcom.almin = alj
        end
        j += 1
    end

    if sd == zero
        if cc > zero
            fn_val = one
            return fn_val
        else
            fn_val = zero
            return fn_val
        end
    end

    sd = sqrt(sd)
    almx = qfcom.almax < -qfcom.almin ? -qfcom.almin : qfcom.almax

    utx = 16 / sd
    up = 4.5 / sd
    un = -up

    findu(n, alb, anc, Ref(utx), 0.5 * acc1)

    if cc != zero && almx > 0.07 * sd
        tausq = 0.25 * acc1 / cfe(n, alb, anc, ith, cc)
        if qfcom.fail
            qfcom.fail = false
        else
            if truncn(n, alb, anc, utx, tausq) < 0.2 * acc1
                qfcom.sigsq += tausq
                findu(n, alb, anc, Ref(utx), 0.25 * acc1)
                trace[6] = sqrt(tausq)
            end
        end
    end

    trace[5] = utx
    acc1 *= 0.5

    d1 = ctff(n, alb, anc, acc1, Ref(up)) - cc
    if d1 < zero
        fn_val = one
        return fn_val
    end

    d2 = cc - ctff(n, alb, anc, acc1, Ref(un))
    if d2 < zero
        fn_val = zero
        return fn_val
    end

    aintv = d1 > d2 ? d1 : d2
    aintv = 2 * pi / aintv

    xnt = utx / aintv
    xntm = 3 / sqrt(acc1)

    if xnt > 1.5 * xntm
        if xntm > xlim
            ifault[] = 1
            return fn_val
        end

        ntm = iround(xntm)
        aintv1 = utx / xntm
        x = 2 * pi / aintv1
        if x <= abs(cc)
            tausq = cfe(n, alb, anc, ith, cc - x) + cfe(n, alb, anc, ith, cc + x)
            tausq = 0.33 * acc1 / (1.1 * tausq)
            if qfcom.fail
                qfcom.fail = false
            else
                acc1 *= 0.67
            end

            integr(n, alb, anc, ntm, aintv1, tausq, false)
            xlim -= xntm
            qfcom.sigsq += tausq
            trace[3] += 1
            trace[2] += ntm + 1

            findu(n, alb, anc, Ref(utx), 0.25 * acc1)
            acc1 *= 0.75
        end
    end

    trace[4] = aintv
    if xnt > xlim
        ifault[] = 1
        return fn_val
    end

    nt = iround(xnt)
    integr(n, alb, anc, nt, aintv, zero, true)
    trace[3] += 1
    trace[2] += nt + 1
    fn_val = 0.5 - qfcom.aintl
    trace[1] = qfcom.ersm
    up = qfcom.ersm

    x = up + acc / 10
    j = 1
    while j <= 8
        if j * x == j * up
            ifault[] = 2
        end
        j *= 2
    end

    trace[7] = qfcom.icount
    return fn_val
end

function countr()
    qfcom.icount += 1
    if qfcom.icount > qfcom.lim
        println("qf: cannot locate integration parameters")
        error("Integration parameters error")
    end
end

function alog1(x::Float64, first::Bool)
    pt1 = 0.1
    one = 1.0
    two = 2.0
    three = 3.0

    if abs(x) > pt1
        if first
            return log(one + x)
        else
            return log(one + x) - x
        end
    end

    y = x / (two + x)
    term = two * y^3
    ak = three
    if first
        s = two
    else
        s = -x
    end

    s = s * y
    y = y^2
    s1 = s + term / ak

    while s1 != s
        ak += two
        term *= y
        s = s1
        s1 = s + term / ak
    end

    return s
end

function exp1(x::Float64)
    neg50 = -50.0
    if x < neg50
        return zero
    else
        return exp(x)
    end
end

function order(alb::Vector{Float64}, ith::Vector{Int})
    for j in 1:qfcom.ir
        alj = abs(alb[j])
        k = j - 1
        while k > 0 && alj > abs(alb[ith[k]])
            ith[k + 1] = ith[k]
            k -= 1
        end
        ith[k + 1] = j
    end
    qfcom.ndtsrt = false
end

function errbd(n::Vector{Int}, alb::Vector{Float64}, anc::Vector{Float64}, uu::Float64, cx::Ref{Float64})
    u = uu
    const = u * qfcom.sigsq
    sum1 = u * const
    u *= two
    for j in 1:qfcom.ir
        nj = n[j]
        alj = alb[j]
        ancj = anc[j]
        x = u * alj
        y = one - x
        const += alj * (ancj / y + nj) / y
        sum1 += ancj * (x / y)^2
        sum1 += nj * (x^2 / y + alog1(-x, false))
    end
    fn_val = exp1(-0.5 * sum1)
    cx[] = const
    return fn_val
end

function ctff(n::Vector{Int}, alb::Vector{Float64}, anc::Vector{Float64}, accx::Float64, upn::Ref{Float64})
    u2 = upn[]
    u1 = zero
    c1 = qfcom.amean

    if u2 > zero
        rb = qfcom.almax
    else
        rb = qfcom.almin
    end

    rb *= two
    u = u2 / (one + u2 * rb)

    while errbd(n, alb, anc, u, Ref{Float64}(c2)) > accx
        u1 = u2
        c1 = c2
        u2 *= two
        u = u2 / (one + u2 * rb)
    end

    u = (c1 - qfcom.amean) / (c2 - qfcom.amean)
    while u < 0.9
        u = (u1 + u2) / two
        if errbd(n, alb, anc, u / (one + u * rb), Ref{Float64}(const)) > accx
            u1 = u
            c1 = const
        else
            u2 = u
            c2 = const
        end
        u = (c1 - qfcom.amean) / (c2 - qfcom.amean)
    end

    fn_val = c2
    upn[] = u2
    return fn_val
end

function truncn(n::Vector{Int}, alb::Vector{Float64}, anc::Vector{Float64}, uu::Float64, tausq::Float64)
    u = uu
    sum1 = zero
    prod2 = zero
    prod3 = zero
    ns = 0
    sum2 = (qfcom.sigsq + tausq) * u^2
    prod1 = two * sum2
    u *= two

    for j in 1:qfcom.ir
        alj = alb[j]
        ancj = anc[j]
        nj = n[j]
        x = (u * alj)^2
        sum1 += ancj * x / (one + x)
        if x > one
            prod2 += nj * log(x)
            prod3 += nj * alog1(x, true)
            ns += nj
        else
            prod1 += nj * alog1(x, true)
        end
    end

    sum1 *= 0.5
    prod2 = prod1 + prod2
    prod3 = prod1 + prod3
    x = exp1(-sum1 - 0.25 * prod2) / pi
    y = exp1(-sum1 - 0.25 * prod3) / pi
    if ns == 0
        err1 = one
    else
        err1 = x * two / ns
    end

    if prod3 > one
        err2 = 2.5 * y
    else
        err2 = one
    end

    if err2 < err1
        err1 = err2
    end

    x = 0.5 * sum2
    if x <= y
        err2 = one
    else
        err2 = y / x
    end

    if err1 < err2
        fn_val = err1
    else
        fn_val = err2
    end

    return fn_val
end

function findu(n::Vector{Int}, alb::Vector{Float64}, anc::Vector{Float64}, utx::Ref{Float64}, accx::Float64)
    ut = utx[]
    u = ut / four
    if truncn(n, alb, anc, u, zero) > accx
        u = ut
        while truncn(n, alb, anc, u, zero) > accx
            ut *= four
            u = ut
        end
    else
        ut = u
        u /= four
        while truncn(n, alb, anc, u, zero) <= accx
            ut = u
            u /= four
        end
    end

    for i in 1:4
        u = ut / divis[i]
        if truncn(n, alb, anc, u, zero) <= accx
            ut = u
        end
    end
    utx[] = ut
end

function integr(n::Vector{Int}, alb::Vector{Float64}, anc::Vector{Float64}, nterm::Int, aintrv::Float64, tausq::Float64, main::Bool)
    ainpi = aintrv / pi
    for k in 0:nterm
        u = (k + 0.5) * aintrv
        sum1 = -two * u * qfcom.c
        sum2 = abs(sum1)
        sum3 = -0.5 * qfcom.sigsq * u^2

        for j in 1:qfcom.ir
            nj = n[j]
            x = two * alb[j] * u
            y = x^2
            sum3 -= 0.25 * nj * alog1(y, true)
            y = anc[j] * x / (one + y)
            z = nj * atan(x) + y
            sum1 += z
            sum2 += abs(z)
            sum3 -= 0.5 * x * y
        end

        x = ainpi * exp1(sum3) / u
        if !main
            x *= one - exp1(-0.5 * tausq * u^2)
        end
        sum1 = sin(0.5 * sum1) * x
        sum2 = 0.5 * sum2 * x
        qfcom.aintl += sum1
        qfcom.ersm += sum2
    end
end

function cfe(n::Vector{Int}, alb::Vector{Float64}, anc::Vector{Float64}, ith::Vector{Int}, x::Float64)
    axl = abs(x)
    sxl = sign(x)
    sum1 = zero
    j = qfcom.ir

    if qfcom.ndtsrt
        order(alb, ith)
    end

    while j > 0
        it = ith[j]
        if alb[it] * sxl > zero
            alj = abs(alb[it])
            axl1 = axl - alj * (n[it] + anc[it])
            axl2 = alj / aln28
            if axl1 > axl2
                axl = axl1
            else
                if axl > axl2
                    axl = axl2
                end
                sum1 = (axl - axl1) / alj
                for k in 1:j-1
                    itk = ith[k]
                    sum1 += n[itk] + anc[itk]
                end
                break
            end
        end
        j -= 1
    end

    if sum1 > 100.0
        fn_val = one
        qfcom.fail = true
    else
        fn_val = 2^(sum1 / four) / (pi * axl^2)
    end

    return fn_val
end

end
