module StudentizedRange

using Printf

# Function to evaluate the tail area of the standardised normal curve
function alnorm(x::Float64, upper::Bool)::Float64
    zero = 0.0
    one = 1.0
    half = 0.5
    con = 1.28

    ltone = 7.0
    utzero = 18.66
    p = 0.398942280444
    q = 0.39990348504
    r = 0.398942280385
    a1 = 5.75885480458
    a2 = 2.62433121679
    a3 = 5.92885724438
    b1 = -29.8213557807
    b2 = 48.6959930692
    c1 = -3.8052e-8
    c2 = 3.98064794e-4
    c3 = -0.151679116635
    c4 = 4.8385912808
    c5 = 0.742380924027
    c6 = 3.99019417011
    d1 = 1.00000615302
    d2 = 1.98615381364
    d3 = 5.29330324926
    d4 = -15.1508972451
    d5 = 30.789933034

    up = upper
    z = x
    if z < zero
        up = !up
        z = -z
    end

    if z <= ltone || (up && z <= utzero)
        y = half * z * z
        if z > con
            fn_val = r * exp(-y) / (z + c1 + d1 / (z + c2 + d2 / (z + c3 + d3 / (z + c4 + d4 / (z + c5 + d5 / (z + c6))))))
        else
            fn_val = half - z * (p - q * y / (y + a1 + b1 / (y + a2 + b2 / (y + a3))))
        end
    else
        fn_val = zero
    end

    if !up
        fn_val = one - fn_val
    end

    return fn_val
end

# Subroutine to produce the normal deviate Z corresponding to a given lower tail area of P
function ppnd16(p::Float64)
    zero = 0.0
    one = 1.0
    half = 0.5
    split1 = 0.425
    split2 = 5.0
    const1 = 0.180625
    const2 = 1.6
    q, r = 0.0, 0.0

    # Coefficients for P close to 0.5
    a0, a1, a2, a3, a4, a5, a6, a7 = 3.3871328727963666080, 1.3314166789178437745e2, 1.9715909503065514427e3, 1.3731693765509461125e4, 4.5921953931549871457e4, 6.7265770927008700853e4, 3.3430575583588128105e4, 2.5090809287301226727e3
    b1, b2, b3, b4, b5, b6, b7 = 4.2313330701600911252e1, 6.8718700749205790830e2, 5.3941960214247511077e3, 2.1213794301586595867e4, 3.9307895800092710610e4, 2.8729085735721942674e4, 5.2264952788528545610e3

    # Coefficients for P not close to 0, 0.5 or 1
    c0, c1, c2, c3, c4, c5, c6, c7 = 1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550, 3.64784832476320460504, 1.27045825245236838258, 2.41780725177450611770e-1, 2.27238449892691845833e-2, 7.74545014278341407640e-4
    d1, d2, d3, d4, d5, d6, d7 = 2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1, 1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9

    # Coefficients for P near 0 or 1
    e0, e1, e2, e3, e4, e5, e6, e7 = 6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580, 2.96560571828504891230e-1, 2.65321895265761230930e-2, 1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7
    f1, f2, f3, f4, f5, f6, f7 = 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2, 7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15

    ifault = 0
    q = p - half
    if abs(q) <= split1
        r = const1 - q * q
        normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
        return normal_dev, ifault
    else
        if q < zero
            r = p
        else
            r = 1.0 - p
        end

        if r <= zero
            ifault = 1
            normal_dev = zero
            return normal_dev, ifault
        end

        r = sqrt(-log(r))
        if r <= split2
            r = r - const2
            normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
        else
            r = r - split2
            normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
        end

        if q < zero
            normal_dev = -normal_dev
        end

        return normal_dev, ifault
    end
end

# Function to calculate an initial quantile for a studentized range distribution
function qtrng0(p::Float64, v::Float64, r::Float64)::Float64
    vmax = 120.0
    half = 0.5
    one = 1.0
    four = 4.0
    c1 = 0.8843
    c2 = 0.2368
    c3 = 1.214
    c4 = 1.208
    c5 = 1.4142

    t, ifault = ppnd16(half + half * p)
    if v < vmax
        t = t + (t^3 + t) / v / four
    end

    q = c1 - c2 * t
    if v < vmax
        q = q - c3 / v + c4 * t / v
    end

    return t * (q * log(r - one) + c5)
end

# Function to evaluate the probability from 0 to q for a studentized range
function prtrng(q::Float64, v::Float64, r::Float64)::Float64
    vw = Vector{Float64}(undef, 30)
    qw = Vector{Float64}(undef, 30)
    zero = 0.0
    fifth = 0.2
    half = 0.5
    one = 1.0
    two = 2.0
    step = 0.45
    pcutj = 0.00003
    pcutk = 0.0001
    vmax = 120.0
    cv1 = 0.193064705
    cv2 = 0.293525326
    cvmax = 0.39894228
    cv = [0.318309886, -0.268132716e-2, 0.347222222e-2, 0.833333333e-1]
    jmin = 3
    jmax = 15
    kmin = 7
    kmax = 15

    prob = zero
    ifault = 0

    if v < one || r < two
        ifault = 1
    end

    if q <= zero || ifault == 1
        return zero
    end

    g = step * r^(-fifth)
    gmid = half * log(r)
    r1 = r - one
    c = log(r * g * cvmax)

    if v <= vmax
        h = step * v^(-half)
        v2 = v * half

        if v == one
            c = cv1
        elseif v == two
            c = cv2
        else
            c = sqrt(v2) * cv[1] / (one + ((cv[2] / v2 + cv[3]) / v2 + cv[4]) / v2)
        end

        c = log(c * r * g * h)
    end

    gstep = g
    qw[1] = -one
    qw[jmax + 1] = -one
    pk1 = one
    pk2 = one

    for k in 1:kmax
        gstep = gstep - g
        gstep = -gstep
        gk = gmid + gstep
        pk = zero

        if pk2 <= pcutk && k > kmin
            continue
        end

        w0 = c - gk * gk * half
        pz = alnorm(gk, true)
        x = alnorm(gk - q, true) - pz

        if x > zero
            pk = exp(w0 + r1 * log(x))
        end

        if v <= vmax
            jump = -jmax

            while true
                jump = jump + jmax

                for j in 1:jmax
                    jj = j + jump

                    if qw[jj] > zero
                        continue
                    end

                    hj = h * j

                    if j < jmax
                        qw[jj + 1] = -one
                    end

                    ehj = exp(hj)
                    qw[jj] = q * ehj
                    vw[jj] = v * (hj + half - ehj * ehj * half)

                    pj = zero
                    x = alnorm(gk - qw[jj], true) - pz

                    if x > zero
                        pj = exp(w0 + vw[jj] + r1 * log(x))
                    end

                    pk = pk + pj

                    if pj > pcutj
                        continue
                    end

                    if jj > jmin || k > kmin
                        break
                    end
                end

                h = -h

                if h >= zero
                    break
                end
            end
        end

        prob = prob + pk

        if k > kmin && pk <= pcutk && pk1 <= pcutk
            break
        end

        pk2 = pk1
        pk1 = pk

        if gstep > zero
            continue
        end
    end

    return prob
end

# Function to approximate the quantile for a studentized range distribution
function qtrng(p::Float64, v::Float64, r::Float64)::Float64
    jmax = 8
    pcut = 0.001
    p75 = 0.75
    p80 = 0.80
    p90 = 0.9
    p99 = 0.99
    p995 = 0.995
    p175 = 1.75
    one = 1.0
    two = 2.0
    five = 5.0
    eps = 1.0e-4

    ifault = 0
    nfault = 0

    if v < one || r < two
        ifault = 1
    end

    if p < p90 || p > p99
        ifault = 2
    end

    if ifault != 0
        return 0.0
    end

    q1 = qtrng0(p, v, r)
    p1 = prtrng(q1, v, r)

    if nfault != 0
        return 0.0
    end

    quantile = q1

    if abs(p1 - p) < pcut
        return quantile
    end

    if p1 > p
        p1 = p175 * p - p75 * p1
    end

    if p1 < p
        p2 = p + (p - p1) * (one - p) / (one - p1) * p75
    end

    if p2 < p80
        p2 = p80
    end

    if p2 > p995
        p2 = p995
    end

    q2 = qtrng0(p2, v, r)

    if nfault != 0
        return 0.0
    end

    for j in 2:jmax
        p2 = prtrng(q2, v, r)

        if nfault != 0
            return 0.0
        end

        e1 = p1 - p
        e2 = p2 - p
        quantile = (q1 + q2) / two
        d = e2 - e1

        if abs(d) > eps
            quantile = (e2 * q1 - e1 * q2) / d
        end

        if abs(e1) < abs(e2)
            q1 = q2
            p1 = p2
        end

        if abs(p1 - p) < pcut * five
            return quantile
        end

        q2 = quantile
    end

    if nfault != 0
        ifault = 9
    end

    return quantile
end

end # module StudentizedRange

# Example usage
function main()
    q = 2.0
    v = 10.0
    r = 5.0

    p = StudentizedRange.prtrng(q, v, r)
    println("Probability from 0 to q: ", p)

    p_quantile = 0.95
    quantile = StudentizedRange.qtrng(p_quantile, v, r)
    println("Quantile for p: ", quantile)
end

main()
