module Sqmcor

using SpecialFunctions  # For loggamma and betainc

# Function to compute the log-gamma function (algorithm AS245)
function alngam(xvalue::Float64)::Float64
    # Coefficients of rational functions
    r1 = [-2.66685511495, -24.4387534237, -21.9698958928, 11.1667541262, 3.13060547623, 0.607771387771, 11.9400905721, 31.4690115749, 15.2346874070]
    r2 = [-78.3359299449, -142.046296688, 137.519416416, 78.6994924154, 4.16438922228, 47.0668766060, 313.399215894, 263.505074721, 43.3400022514]
    r3 = [-212159.572323, 230661.510616, 27464.764705, -40262.1119975, -2296.60729780, -116328.495004, -146025.937511, -24235.7409629, -570.691009324]
    r4 = [0.279195317918525, 0.4917317610505968, 0.0692910599291889, 3.350343815022304, 6.012459259764103]

    # Fixed constants
    alr2pi = 0.918938533204673
    four = 4.0
    half = 0.5
    one = 1.0
    onep5 = 1.5
    twelve = 12.0
    zero = 0.0

    # Machine-dependent constants
    xlge = 5.10e6
    xlgst = floatmax(Float64)

    x = xvalue
    fn_val = zero

    # Test for valid function argument
    if x >= xlgst
        println("AS 245: Argument x too large")
        return fn_val
    end

    if x <= zero
        println("AS 245: Argument x <= 0")
        return fn_val
    end

    # Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined
    if x < onep5
        if x < half
            fn_val = -log(x)
            y = x + one

            # Test whether X < machine epsilon
            if y == one
                return fn_val
            end
        else
            fn_val = zero
            y = x
            x = (x - half) - half
        end
        fn_val += x * (((r1[5]*y + r1[4])*y + r1[3])*y + r1[2])*y + r1[1] / (((y + r1[9])*y + r1[8])*y + r1[7])*y + r1[6]
        return fn_val
    end

    # Calculation for 1.5 <= X < 4.0
    if x < four
        y = (x - one) - one
        fn_val = y * (((r2[5]*x + r2[4])*x + r2[3])*x + r2[2])*x + r2[1] / (((x + r2[9])*x + r2[8])*x + r2[7])*x + r2[6]
        return fn_val
    end

    # Calculation for 4.0 <= X < 12.0
    if x < twelve
        fn_val = (((r3[5]*x + r3[4])*x + r3[3])*x + r3[2])*x + r3[1] / (((x + r3[9])*x + r3[8])*x + r3[7])*x + r3[6]
        return fn_val
    end

    # Calculation for X >= 12.0
    y = log(x)
    fn_val = x * (y - one) - half * y + alr2pi
    if x > xlge
        return fn_val
    end
    x1 = one / x
    x2 = x1 * x1
    fn_val += x1 * ((r4[3]*x2 + r4[2])*x2 + r4[1]) / ((x2 + r4[5])*x2 + r4[4])
    return fn_val
end

# Function to compute the incomplete beta function ratio (algorithm AS63)
function betain(x::Float64, p::Float64, q::Float64, beta::Float64)::Float64
    # Local variables
    indx = false
    psq = p + q
    cx = 1.0 - x
    acu = 1.0e-14

    if p <= 0.0 || q <= 0.0
        println("AS63: Either p or q <= 0")
        return x
    end

    if x < 0.0 || x > 1.0
        println("AS63: Argument x outside range (0, 1)")
        return x
    end

    if x == 0.0 || x == 1.0
        return x
    end

    if p < psq * x
        xx = cx
        cx = x
        pp = q
        qq = p
        indx = true
    else
        xx = x
        pp = p
        qq = q
        indx = false
    end

    term = 1.0
    ai = 1.0
    fn_val = 1.0
    ns = floor(Int, qq + cx * psq)

    # Use Soper's reduction formulae.
    rx = xx / cx
    while true
        temp = qq - ai
        if ns == 0
            rx = xx
        end
        term *= temp * rx / (pp + ai)
        fn_val += term
        temp = abs(term)
        if temp <= acu && temp <= acu * fn_val
            break
        end
        ai += 1.0
        ns -= 1
        if ns < 0
            temp = psq
            psq += 1.0
        end
    end

    fn_val *= exp(pp * log(xx) + (qq - 1.0) * log(cx) - beta) / pp
    if indx
        fn_val = 1.0 - fn_val
    end

    return fn_val
end

# Function to compute the CDF for the distribution of the square of the multiple correlation coefficient (algorithm AS260)
function sqmcor(x::Float64, ip::Int, n::Int, rho2::Float64)::Float64
    errmax = 1.0e-8
    zero = 0.0
    half = 0.5
    one = 1.0
    itrmax = 100

    fn_val = x
    if rho2 < zero || rho2 > one || ip < 2 || n <= ip
        println("Error in Function sqmcor, rho2 must be in the range (0,1)")
        return fn_val
    end

    if x < zero || x > one
        println("Error in Function sqmcor, x must be in the range (0,1)")
        return fn_val
    end

    if x == zero || x == one
        return fn_val
    end

    a = half * (ip - 1)
    b = half * (n - ip)

    # Initialize the series
    beta = exp(alngam(a) + alngam(b) - alngam(a + b))
    temp = betain(x, a, b, beta)

    gx = exp(a * log(x) + b * log(one - x) - log(a)) / beta
    q = (one - rho2) ^ (a + b)
    xj = zero
    term = q * temp
    sumq = one - q
    fn_val = term

    # Perform recurrence until convergence is achieved
    while true
        xj += one
        temp -= gx
        gx = gx * (a + b + xj - one) * x / (a + xj)
        q = q * (a + b + xj - one) * rho2 / xj
        sumq -= q
        term = temp * q
        fn_val += term

        # Check for convergence and act accordingly
        errbd = (temp - gx) * sumq
        if xj < itrmax && errbd > errmax
            continue
        end

        if errbd > errmax
            println("Error in SQMCOR, errbd > errmax")
        end
        break
    end

    return fn_val
end

end  # module Sqmcor

# Example usage of the sqmcor function
using .Sqmcor

function main()
    x = 0.5
    ip = 3
    n = 30
    rho2 = 0.1
    result = Sqmcor.sqmcor(x, ip, n, rho2)
    println("The CDF value for the given parameters is $result")
end

main()
