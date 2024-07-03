module ShapiroWilk

using Printf, Random, Statistics

# The main function that calculates the Shapiro-Wilk W test and its significance level
function swilk!(init::Bool, x::Vector{Float64}, n::Int, n1::Int, n2::Int, a::Vector{Float64})
    w = 0.0
    pw = 1.0
    ifault = 0

    upper = true
    c1 = [0.0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056]
    c2 = [0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633]
    c3 = [0.5440, -0.39978, 0.025054, -0.6714e-3]
    c4 = [1.3822, -0.77857, 0.062767, -0.0020322]
    c5 = [-1.5861, -0.31082, -0.083751, 0.0038915]
    c6 = [-0.4803, -0.082676, 0.0030302]
    c7 = [0.164, 0.533]
    c8 = [0.1736, 0.315]
    c9 = [0.256, -0.00635]
    g = [-2.273, 0.459]
    z90 = 1.2816
    z95 = 1.6449
    z99 = 2.3263
    zm = 1.7509
    zss = 0.56268
    bf1 = 0.8378
    xx90 = 0.556
    xx95 = 0.622
    zero = 0.0
    one = 1.0
    two = 2.0
    three = 3.0
    sqrth = 0.70711
    qtr = 0.25
    th = 0.375
    small = 1e-19
    pi6 = 1.909859
    stqr = 1.047198

    pw = one
    if w >= zero
        w = one
    end

    an = n
    ifault = 3
    nn2 = div(n, 2)
    if n2 < nn2
        return w, pw, ifault
    end

    ifault = 1
    if n < 3
        return w, pw, ifault
    end

    # If INIT is false, calculates coefficients for the test
    if !init
        if n == 3
            a[1] = sqrth
        else
            an25 = an + qtr
            summ2 = zero
            for i in 1:n2
                ai, ifault = ppnd7((i - th) / an25)
                a[i] = ai
                summ2 += a[i] ^ 2
            end
            summ2 *= two
            ssumm2 = sqrt(summ2)
            rsn = one / sqrt(an)
            a1 = poly(c1, rsn) - a[1] / ssumm2
            
            # Normalize coefficients
            if n > 5
                i1 = 3
                a2 = -a[2] / ssumm2 + poly(c2, rsn)
                fac = sqrt((summ2 - two * a[1] ^ 2 - two * a[2] ^ 2) / (one - two * a1 ^ 2 - two * a2 ^ 2))
                a[1] = a1
                a[2] = a2
            else
                i1 = 2
                fac = sqrt((summ2 - two * a[1] ^ 2) / (one - two * a1 ^ 2))
                a[1] = a1
            end
            for i in i1:nn2
                a[i] = -a[i] / fac
            end
        end
        init = true
    end

    if n1 < 3
        return w, pw, ifault
    end

    ncens = n - n1
    ifault = 4
    if ncens < 0 || (ncens > 0 && n < 20)
        return w, pw, ifault
    end

    ifault = 5
    delta = float(ncens) / an
    if delta > 0.8
        return w, pw, ifault
    end

    # If W input as negative, calculate significance level of -W
    if w < zero
        w1 = one + w
        ifault = 0
        return one - w1, alnorm(log(w1) / poly(c3, log(an)), upper), ifault
    end

    # Check for zero range
    ifault = 6
    range = x[n1] - x[1]
    if range < small
        return w, pw, ifault
    end

    # Check for correct sort order on range - scaled X
    ifault = 7
    xx = x[1] / range
    sx = xx
    sa = -a[1]
    j = n - 1
    for i in 2:n1
        xi = x[i] / range
        if xx - xi > small
            println("x(i)s out of order")
            return w, pw, ifault
        end
        sx += xi
        if i != j
            sa += sign(1, i - j) * a[min(i, j)]
        end
        xx = xi
        j -= 1
    end
    ifault = 0
    if n > 5000
        ifault = 2
    end

    # Calculate W statistic as squared correlation between data and coefficients
    sa /= n1
    sx /= n1
    ssa = zero
    ssx = zero
    sax = zero
    j = n
    for i in 1:n1
        if i != j
            asa = sign(1, i - j) * a[min(i, j)] - sa
        else
            asa = -sa
        end
        xsx = x[i] / range - sx
        ssa += asa ^ 2
        ssx += xsx ^ 2
        sax += asa * xsx
        j -= 1
    end

    # W1 equals (1-W) calculated to avoid excessive rounding error for W very near 1
    ssassx = sqrt(ssa * ssx)
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx)
    w = one - w1

    # Calculate significance level for W (exact for N=3)
    if n == 3
        pw = pi6 * (asin(sqrt(w)) - stqr)
        return w, pw, ifault
    end

    y = log(w1)
    xx = log(an)
    m = zero
    s = one
    if n <= 11
        gamma = poly(g, an)
        if y >= gamma
            pw = small
            return w, pw, ifault
        end
        y = -log(gamma - y)
        m = poly(c3, an)
        s = exp(poly(c4, an))
    else
        m = poly(c5, xx)
        s = exp(poly(c6, xx))
    end

    if ncens > 0
        # Censoring by proportion NCENS/N. Calculate mean and sd of normal equivalent deviate of W.
        ld = -log(delta)
        bf = one + xx * bf1
        z90f = z90 + bf * poly(c7, xx90 ^ xx) * ld
        z95f = z95 + bf * poly(c8, xx95 ^ xx) * ld
        z99f = z99 + bf * poly(c9, xx) * ld

        # Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get pseudo-mean and pseudo-sd of z as the slope and intercept
        zfm = (z90f + z95f + z99f) / three
        zsd = (z90 * (z90f - zfm) + z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss
        zbar = zfm - zsd * zm
        m += zbar * s
        s *= zsd
    end

    pw = alnorm((y - m) / s, upper)
    return w, pw, ifault
end

# Subroutine to produce the normal deviate Z corresponding to a given lower tail area of P
function ppnd7(p::Float64)
    low_prec = Float32
    zero = 0.0
    one = 1.0
    half = 0.5
    split1 = 0.425
    split2 = 5.0
    const1 = 0.180625
    const2 = 1.6
    q, r = 0.0, 0.0

    # Coefficients for P close to 0.5
    a0, a1, a2, a3 = 3.3871327179e0, 5.0434271938e1, 1.5929113202e2, 5.9109374720e1
    b1, b2, b3 = 1.7895169469e1, 7.8757757664e1, 6.7187563600e1

    # Coefficients for P not close to 0, 0.5 or 1
    c0, c1, c2, c3 = 1.4234372777e0, 2.7568153900e0, 1.3067284816e0, 1.7023821103e-1
    d1, d2 = 7.3700164250e-1, 1.2021132975e-1

    # Coefficients for P near 0 or 1
    e0, e1, e2, e3 = 6.6579051150e0, 3.0812263860e0, 4.2868294337e-1, 1.7337203997e-2
    f1, f2 = 2.4197894225e-1, 1.2258202635e-2

    ifault = 0
    q = p - half
    if abs(q) <= split1
        r = const1 - q * q
        normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / (((b3 * r + b2) * r + b1) * r + one)
        return normal_dev, ifault
    else
        if q < zero
            r = p
        else
            r = one - p
        end

        if r <= zero
            ifault = 1
            normal_dev = zero
            return normal_dev, ifault
        end

        r = sqrt(-log(r))
        if r <= split2
            r = r - const2
            normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one)
        else
            r = r - split2
            normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one)
        end

        if q < zero
            normal_dev = -normal_dev
        end

        return normal_dev, ifault
    end
end

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

# Function to calculate the algebraic polynomial of order nord-1 with array of coefficients c
function poly(c::Vector{Float64}, x::Float64)::Float64
    nord = length(c)
    fn_val = c[1]
    if nord == 1
        return fn_val
    end

    p = x * c[nord]
    if nord == 2
        return fn_val + p
    end

    n2 = nord - 2
    j = n2 + 1
    for i in 1:n2
        p = (p + c[j]) * x
        j -= 1
    end

    return fn_val + p
end

# Function to generate normally distributed pseudo-random number with zero mean and unit variance
function random_normal()
    s = 0.449871
    t = -0.386595
    a = 0.19600
    b = 0.25472
    half = 0.5
    r1 = 0.27597
    r2 = 0.27846

    while true
        u = rand()  # Uniform random number between 0 and 1
        v = 1.7156 * (rand() - half)
        
        x = u - s
        y = abs(v) - t
        q = x^2 + y * (a * y - b * x)
        
        if q < r1
            break
        end
        
        if q > r2
            continue
        end
        
        if v^2 < -4.0 * log(u) * u^2
            break
        end
    end

    return v / u
end

# Function to sort list and return the order of elements
function quick_sort!(list::Vector{Float64})
    order = collect(1:length(list))
    quick_sort!(list, order, 1, length(list))
    return order
end

# Recursive helper function for quick_sort!
function quick_sort!(list::Vector{Float64}, order::Vector{Int}, left::Int, right::Int)
    if left < right
        pivotIndex = partition!(list, order, left, right)
        quick_sort!(list, order, left, pivotIndex - 1)
        quick_sort!(list, order, pivotIndex + 1, right)
    end
end

# Partition function for quick_sort!
function partition!(list::Vector{Float64}, order::Vector{Int}, left::Int, right::Int)::Int
    pivot = list[right]
    i = left - 1
    for j in left:right-1
        if list[j] <= pivot
            i += 1
            list[i], list[j] = list[j], list[i]
            order[i], order[j] = order[j], order[i]
        end
    end
    list[i+1], list[right] = list[right], list[i+1]
    order[i+1], order[right] = order[right], order[i+1]
    return i + 1
end

end # module ShapiroWilk

# Example usage
function main()
    n = 10
    n1 = n
    n2 = div(n, 2)
    x = [ShapiroWilk.random_normal() for _ in 1:n1]
    a = Vector{Float64}(undef, n2)
    init = false

    # Sort the X's into ascending order
    order = ShapiroWilk.quick_sort!(x)

    # Apply the Shapiro-Wilk test
    w, pw, ifault = ShapiroWilk.swilk!(init, x, n, n1, n2, a)

    if ifault != 0
        println("IFault = $ifault")
    end

    println("Shapiro-Wilk W = $w")
    println("P-value = $pw")
end

main()
