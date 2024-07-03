module Rsquared
using Printf

export qr2, ifault

# Module-level variable to track fault status
ifault = 0

# Function to compute the quantile of the distribution of the square of the sample multiple correlation coefficient
function qr2(m::Int, size::Int, rho2::Float64, p::Float64)::Float64
    # Local variables
    zero = 0.0
    half = 0.5
    one = 1.0
    two = 2.0
    eps = 1.0e-06
    delta = 1.0e-04
    rp = 1.772453850905516028
    itrmax = 10

    quantile = p

    # Test for admissibility of arguments
    global ifault
    ifault = 2
    if m <= 1 || size <= m || rho2 < zero || rho2 > one || p < zero || p > one
        return quantile
    end
    ifault = 0
    if p == zero || p == one
        return quantile
    end

    # Calculate the constants needed for each Newton's iteration
    a = (m - 1) / two
    b = (size - m) / two
    ab = (size - 1) / two

    if (m + 1) % 2 == 0
        na = a + half
        ga = one
        for i in 1:na
            ga *= i
        end
    else
        na = a + one
        ga = rp
        for i in 1:na
            ga *= (i - half)
        end
    end

    if (size - m) % 2 == 0
        nb = b - half
        gb = one
        for i in 1:nb
            gb *= i
        end
    else
        nb = b
        gb = rp
        for i in 1:nb
            gb *= (i - half)
        end
    end

    if (size - 1) % 2 == 0
        nab = ab - half
        gab = one
        for i in 1:nab
            gab *= i
        end
    else
        nab = ab
        gab = rp
        for i in 1:nab
            gab *= (i - half)
        end
    end

    q0 = (one - rho2) ^ ab
    coeff = gab / ga / gb

    # Use 0.5 as a starting value for Newton's iterations
    y = half

    # Perform Newton's iterations
    for iter in 1:itrmax
        # Evaluate the first terms of the series for CDF (distribution function) and PDF (density)
        n = one
        yp = one - y
        t = coeff * y ^ a * yp ^ b
        s = a * t / y / yp
        q = q0
        v = q
        cdf = v * t
        pdf = q * s

        # Check if a + n > (a + b + n)y
        while true
            if a + n > (a + b + n) * y
                break
            end

            # Evaluate the next terms of two series and then the partial sums
            q = q * (a + b + n - one) * rho2 / n
            v += q
            s = t * (a + b + n - one) / yp
            t = t * y * (a + b + n - one) / (a + n)
            cdf += v * t
            pdf += q * s
            n += one
        end

        # Find the error bounds and check for convergence for both series
        bndcdf = t * y * (a + b + n - one) / (a + n - (a + b + n) * y)
        bndpdf = t * (a + b + n - one) * (one - v) / yp

        while true
            if bndcdf <= eps && bndpdf <= eps
                break
            end

            # Continue to update the terms and then accumulate
            q = q * (a + b + n - one) * rho2 / n
            v += q
            if bndcdf <= eps
                s = s * y * (a + b + n - one) / (a + n - one)
                pdf += q * s
                n += one
                bndpdf = s * y * (a + b + n - one) * (one - v) / (a + n - one)
            elseif bndpdf <= eps
                t = t * y * (a + b + n - one) / (a + n)
                cdf += v * t
                n += one
                bndcdf = t * y * (a + b + n - one) / (a + n - (a + b + n) * y)
            else
                s = t * (a + b + n - one) / yp
                t = t * y * (a + b + n - one) / (a + n)
                cdf += v * t
                pdf += q * s
                n += one
            end
        end

        # Obtain a new Y and make changes if it is illegal
        diff = (cdf - p) / pdf
        ynew = y - diff
        if ynew <= zero
            y = y / two
        elseif ynew >= one
            y = (y + one) / two
        else
            y = ynew
        end

        # Check for convergence of Newton's iterations
        if abs(diff) <= delta * y
            quantile = y
            return quantile
        end
    end

    ifault = 1
    return quantile
end

end  # module Rsquared

# Main program to call qr2 and produce output
using .Rsquared

function main()
    while true
        println("ENTER M (>1), N (>M), RHO2 (BETWEEN 0 AND 1), and P (BETWEEN 0 AND 1) FOR QR2 ==> ")
        input = readline()
        m, size, rho2, p = split(input)
        m = parse(Int, m)
        size = parse(Int, size)
        rho2 = parse(Float64, rho2)
        p = parse(Float64, p)

        y = Rsquared.qr2(m, size, rho2, p)
        icode = Rsquared.ifault + 1

        if icode == 1
            @printf("qr2(%d, %d, %.4f, %.4f) = %.5f\n", m, size, rho2, p, y)
        elseif icode == 2
            println("NO convergence after 10 Newton ITERATIONS")
        elseif icode == 3
            println("THE INPUT VALUE IS ILLEGAL!")
        end

        println("ENTER 1 TO CONTINUE OR 0 TO QUIT ==> ")
        k = parse(Int, readline())
        if k != 1
            break
        end
    end
end

main()
