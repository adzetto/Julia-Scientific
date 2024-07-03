module Mult_Normal

using LinearAlgebra

# Define global variables
global ix, iy, iz

export mulnor

const dp = Float64

function mulnor(covar::Vector{dp}, n::Int, f::Function, sdev::dp, rhigh::dp, init::Int, iseed::Int, t::Vector{dp}, work::Matrix{dp})
    """
    Finds the probability that a normally distributed random N-vector with mean 0 and covariance COVAR falls in an area enclosed by the external user-defined function F.
    """
    const maxitr = 10000
    const zero = 0.0
    const cox = 1.0 + 2.0 / max(25, min(init, 1000))
    
    prob = zero
    var = zero
    
    # Find Cholesky decomposition of covariance matrix
    nn = n * (n + 1) // 2
    t, nullty, ifault = chol(covar, n, nn, t)
    
    q = n - nullty
    ifault = ifault * ifault
    if ifault != 0
        return prob, ifault
    end
    
    # Transpose Cholesky factor, omitting columns of zeroes. Resulting vector is packed form of T'J, written row-wise beginning with last row.
    work[1:n, 1:n] .= zero
    ii = 0
    for i in 1:n
        for j in 1:i
            ii += 1
            work[i, j] = t[ii]
        end
    end
    
    j = 0
    for i in n:-1:1
        if j == nullty
            break
        end
        if work[i, i] == zero
            j += 1
            for k in i + 1:n
                for kk in i:n - j
                    work[k, kk] = work[k, kk + 1]
                end
            end
        end
    end
    
    k = 0
    for j in n:-1:1
        for i in 1:min(j, q)
            k += 1
            t[k] = work[j, i]
        end
    end
    
    # Take pilot sample; estimate variance of prob from ORTHP.
    for k in 1:init1
        p, iseed, ifault = orthp(f, rhigh, n, q, t, work, iseed)
        if ifault != 0
            return prob, ifault
        end
        prob += p
        var += p * p
    end
    
    prob /= init1
    var = var / init1 - prob * prob
    
    iter = Int(cox * var / (sdev * sdev)) + 1
    if iter > maxitr
        ifault = 5
        return prob, ifault
    end
    
    for k in init1 + 1:iter
        p, iseed, ifault = orthp(f, rhigh, n, q, t, work, iseed)
        if ifault != 0
            return prob, ifault
        end
        prob += (p - prob) / k
    end
    
    return prob, ifault
end

function chol(a::Vector{dp}, n::Int, nn::Int, u::Vector{dp})
    const eta = 1.0e-05
    const zero = 0.0
    
    ifault = 1
    if n <= 0
        return u, 0, ifault
    end
    
    ifault = 3
    if nn != n * (n + 1) // 2
        return u, 0, ifault
    end
    
    ifault = 2
    nullty = 0
    j = 1
    k = 0
    eta2 = eta * eta
    ii = 0
    
    for icol in 1:n
        ii += icol
        x = eta2 * a[ii]
        l = 0
        kk = 0
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
            if u[l] == zero
                u[k] = w / u[l]
                continue
            end
            if w * w > abs(x * a[kk])
                return u, 0, ifault
            end
            u[k] = zero
        end
        
        if abs(w) <= abs(eta * a[k])
            u[k] = zero
            nullty += 1
            continue
        end
        
        if w < zero
            return u, 0, ifault
        end
        
        u[k] = sqrt(w)
    end
    
    ifault = 0
    return u, nullty, ifault
end

function chiprb(x::dp, df::Int)
    const rtpirp = 0.564189583547756
    const big = 88.03
    const small = -85.2
    const zero = 0.0
    const one = 1.0
    const half = 0.5
    
    s = zero
    if x > big
        return s
    end
    
    s = one
    if x < small
        return s
    end
    
    a = half * x
    y = exp(-a)
    odd = mod(df, 2)
    if odd == 0
        s = y
        e = one
        z = zero
    else
        sqrta = sqrt(a)
        s = zerfp1(-sqrta)
        e = rtpirp / sqrta
        z = -half
    end
    
    if df < 3
        return s
    end
    
    n2 = df // 2 - 1 + odd
    c = zero
    for i in 1:n2
        z += 1
        e = e * (a / z)
        c += e
    end
    
    s += c * y
    return s
end

function orthp(f::Function, rhigh::dp, n::Int, q::Int, t::Vector{dp}, z::Matrix{dp}, iseed::Int)
    const zero = 0.0
    const rt2rcp = 0.707106781
    
    prob = zero
    qp1 = q + 1
    qp2 = q + 2
    
    # Put randomly oriented orthonormal system of q-vectors in Z.
    z, fab, iseed = rnortm(n, q, qp1, 1, 1, z[:, qp1], iseed, z)
    
    # Replace each column of Z by TRANS(CHOLESKY FACTOR)*J*COLUMN OF Z.
    for i in 1:q
        jj = 0
        for j in n:-1:1
            s = zero
            for k in 1:min(j, q)
                jj += 1
                s += t[jj] * z[k, i]
            end
            z[j, i] = s
        end
    end
    
    # Estimate the probability that N(0, COVAR) is within the boundary given by F.
    for i in 1:q
        a = zero
        b = rhigh
        r, ifault = zerovc(f, a, b, n, z[:, i])
        prob -= chiprb(r * r, q)
    end
    
    for i in 1:q
        a = zero
        b = rhigh
        for k in 1:n
            z[k, qp1] = -z[k, i]
        end
        r, ifault = zerovc(f, a, b, n, z[:, qp1])
        prob -= chiprb(r * r, q)
    end
    
    for i in 1:q-1
        for j in i+1:q
            for jj in -1:1:2
                for ii in -1:1:2
                    for k in 1:n
                        z[k, qp1] = rt2rcp * (z[k, i] * ii + z[k, j] * jj)
                    end
                    a = zero
                    b = rhigh
                    r, ifault = zerovc(f, a, b, n, z[:, qp1])
                    prob -= chiprb(r * r, q)
                end
            end
        end
    end
    
    prob = prob / (2.0 * q * q) + 1.0
    return prob, iseed, ifault
end

function rnortm(n::Int, np1::Int, ibconf::Int, nb::Int, b::Vector{dp}, iseed::Int, a::Matrix{dp})
    const zero = 0.0
    const one = 1.0
    
    fab = one
    a .= zero
    diagm!(a, one)
    
    chisq = zeros(dp, n)
    for l in 1:n
        np1l = np1 - l
        u, chisq[l] = normal(np1l + 1, iseed)
        ul = sqrt(chisq[l])
        if ul == zero
            return a, fab, iseed
        end
        if l == n
            break
        end
        
        w = copy(u[l:np1])
        wl = sign(ul, -u[l])
        wil = one / wl
        wwil = one / (wl - w[l])
        
        for i in 1:n
            t2 = dot(a[i, l:n], w[l:n])
            hv = t2 * wil
            hw = (hv - a[i, l]) * wwil
            a[i, l] = hv
            for j in l + 1:n
                a[i, j] -= hw * w[j]
            end
        end
        
        if u[l] > zero
            continue
        end
        for i in 1:n
            a[i, l] = -a[i, l]
        end
    end
    
    fab = abs(fab)
    return a, fab, iseed
end

function zerovc(f::Function, r0::dp, r1::dp, n::Int, v::Vector{dp})
    const half = 0.5
    const eta = 1.0e-05
    const niter = 500
    
    y = r0 * v
    f0 = f(n, y)
    y = r1 * v
    f1 = f(n, y)
    if sign(half, f1) == sign(half, f0)
        ifault = 3
        return 0.0, ifault
    end
    
    for j in 1:niter
        r = r1 - f1 * (r1 - r0) / (f1 - f0)
        y = r * v
        fr = f(n, y)
        if abs(r - r1) < eta
            return r, 0
        end
        if sign(half, fr) != sign(half, f1)
            r0 = r1
            f0 = f1
        else
            f0 *= half
        end
        r1 = r
        f1 = fr
    end
    
    return r, 6
end

function zerfp1(a::dp)
    const sqrt2 = 1.41421356237309
    return 2.0 * alnorm(sqrt2 * a, false)
end

function normal(np1::Int, iseed::Int)
    chisq = 0.0
    u = zeros(dp, np1)
    for i in 1:np1 - 1
        u[i], iseed = ppnd(random(iseed))
        chisq += u[i] * u[i]
    end
    return u, chisq
end

function ppnd(p::dp)
    const zero = 0.0
    const half = 0.5
    const one = 1.0
    const split = 0.42
    const a0 = 2.50662823884
    const a1 = -18.61500062529
    const a2 = 41.39119773534
    const a3 = -25.44106049637
    const b1 = -8.47351093090
    const b2 = 23.08336743743
    const b3 = -21.06224101826
    const b4 = 3.13082909833
    const c0 = -2.78718931138
    const c1 = -2.29796479134
    const c2 = 4.85014127135
    const c3 = 2.32121276858
    const d1 = 3.54388924762
    const d2 = 1.63706781897
    
    q = p - half
    if abs(q) > split
        r = ifelse(q > zero, 1 - p, p)
        if r <= zero
            return zero, 1
        end
        r = sqrt(-log(r))
        fn_val = ((c3 * r + c2) * r + c1) * r + c0 / ((d2 * r + d1) * r + one)
        return ifelse(q < zero, -fn_val, fn_val), 0
    end
    
    r = q * q
    fn_val = q * (((a3 * r + a2) * r + a1) * r + a0) / (((b4 * r + b3) * r + b2) * r + b1)
    return fn_val, 0
end

function random(iseed::Int)
    ix = 171 * mod(ix, 177) - 2 * (ix ÷ 177)
    iy = 172 * mod(iy, 176) - 35 * (iy ÷ 176)
    iz = 170 * mod(iz, 178) - 63 * (iz ÷ 178)
    
    ix += ifelse(ix < 0, 30269, 0)
    iy += ifelse(iy < 0, 30307, 0)
    iz += ifelse(iz < 0, 30323, 0)
    
    return mod(Float64(ix) / 30269.0 + Float64(iy) / 30307.0 + Float64(iz) / 30323.0, 1.0)
end

function alnorm(x::dp, upper::Bool)
    const ltone = 7.0
    const utzero = 18.66
    const zero = 0.0
    const half = 0.5
    const one = 1.0
    const con = 1.28
    
    z = abs(x)
    y = zero
    up = upper
    if x < zero
        z = -x
        up = !up
    end
    
    if z <= ltone || up && z <= utzero
        y = 0.5 * z * z
        if z <= con
            return 0.5 - z * (0.398942280444 - 0.399903438504 * y / (y + 5.75885480458 - 29.8213557808 / (y + 2.62433121679 + 48.6959930692 / (y + 5.92885724438))))
        end
        return 0.398942280385 * exp(-y) / (z - 3.8052e-08 + 1.00000615302 / (z + 3.98064794e-4 + 1.98615381364 / (z - 0.151679116635 + 5.29330324926 / (z + 4.8385912808 - 15.1508972451 / (z + 0.742380924027 + 30.789933034 / (z + 3.99019417011))))))
    end
    if !up
        return one - y
    end
    return zero
end

end # module Mult_Normal

using .Mult_Normal

function f1(n, v)
    """
    Boundary of region is an n-dimensional cube.
    CUBSZ is the distance from the origin to the side of the cube.
    """
    const cubsz = 3.0
    return maximum(abs(v)) - cubsz
end

function f2(n, v)
    """
    The boundary is an ellipsoid of the form X'AX + BX + C = 0.
    The upper triangle of A is stored column-wise in packed form.
    In this example, A = inverse([1  .3; .3  1]), B = 0, C = -1.
    """
    const a = [1.0989011, -0.32967033, 1.0989011]
    const b = zeros(10)
    const c = -1.0
    fn_val = c
    k = 0
    for i in 1:n
        k += 1
        fn_val += (v[i]^2) * a[k]
        for j in i+1:n
            k += 1
            fn_val += 2.0 * v[i] * v[j] * a[k]
        end
    end
    return fn_val + dot(b[1:n], v)
end

function f3(n, v)
    """
    Rectangle is given by A < V < B for each coordinate.
    This subroutine may also be used to calculate the CDF by
    making the appropriate sides of the rectangle sufficiently large.
    This example is the four-dimensional rectangle defined
    by -infinity (almost) < V(i) < 2 for each i.
    """
    const a = [-20.0, -20.0, -20.0, -20.0]
    const b = [2.0, 2.0, 2.0, 2.0]
    fn_val = -1.0
    for i in 1:n
        temp = ifelse(v[i] > 0.0, v[i] - b[i], a[i] - v[i])
        fn_val = max(fn_val, temp)
    end
    return fn_val
end

function f4(n, v)
    """
    Boundary is sphere of radius sqrt(radsq).
    """
    const radsq = 4.0
    return sum(v.^2) - radsq
end

function f5(n, v)
    """
    A convex polytope with P faces is defined by the inequalities
    A(1,i)*V(1) + A(2,i)*V(2) +...+ A(n,i)*V(n) <= 1
    for i=1,...,P.
    The figure in this demonstration program is the
    regular pentagon with distance from center to point = 1.
    """
    const p = 5
    const a = [0.726542528 -1.0, 1.175570504, 0.381966011, 0.0, 1.236067978,
               -0.726542528, -1.0, -1.175570504, 0.381966011]
    fn_val = -1.0
    for i in 1:n
        fn_val += a[i, 1] * v[i]
    end
    for i in 2:p
        temp = -1.0
        for j in 1:n
            temp += a[j, i] * v[j]
        end
        fn_val = max(fn_val, temp)
    end
    return fn_val
end

function f6(n, v)
    """
    The figure in this demonstration program is the
    regular pentagram with distance from center to point = 2.618.
    """
    const p = 5
    const tang = [-1.376381921, -0.3249196962329, 0.3249196962329, 1.376381921, 0.0]
    const a = [-0.726542528, -1.0, 1.175570504, 0.381966011, 0.726542528,
               -1.0, 0.0, 1.236067978, 1.175570504, 0.381966011]
    x = [abs(v[1]), v[2]]
    index = p
    if x[1] == 0.0
        index = 3
    else
        tangnt = x[2] / x[1]
        for i in 1:p
            if tangnt <= tang[i]
                index = i
                break
            end
        end
    end
    return -1.0 + dot(a[:, index], x)
end

function main()
    n = 2
    sdev = 1.0
    rhigh = 5.0
    init = 25
    iseed = 1234567
    
    covar = [1.0, 0.3, 1.0]
    work = zeros(dp, 50, 52)
    t = zeros(dp, 1275)
    
    functions = [f1, f2, f3, f4, f5, f6]
    
    for i in 1:6
        prob, ifault = Mult_Normal.mulnor(covar, n, functions[i], sdev, rhigh, init, iseed, t, work)
        println("Function f$i")
        println("Estimated Probability: $prob with IFAULT = $ifault")
    end
end

main()
