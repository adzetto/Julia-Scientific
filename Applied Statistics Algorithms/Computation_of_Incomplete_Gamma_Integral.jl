module GammaDistribution

using SpecialFunctions

# Function to compute the Incomplete Gamma Integral
function gammad(x::Float64, p::Float64)::Float64
    # Constants
    zero = 0.0
    one = 1.0
    two = 2.0
    three = 3.0
    nine = 9.0
    tol = 1.0e-14
    xbig = 1.0e8
    plimit = 1000.0
    elimit = -88.0
    oflo = 1.0e37

    fn_val = zero

    # Check that we have valid values for X and P
    if p <= zero || x < zero
        println("AS239: Either p <= 0 or x < 0")
        return fn_val
    end
    if x == zero
        return fn_val
    end

    # Use a normal approximation if P > PLIMIT
    if p > plimit
        pn1 = three * sqrt(p) * ((x / p)^(one / three) + one / (nine * p) - one)
        fn_val = cdf(Normal(), pn1)
        return fn_val
    end

    # If X is extremely large compared to P then set fn_val = 1
    if x > xbig
        fn_val = one
        return fn_val
    end

    if x <= one || x < p
        # Use Pearson's series expansion.
        arg = p * log(x) - x - loggamma(p + one)
        c = one
        fn_val = one
        a = p
        while true
            a += one
            c *= x / a
            fn_val += c
            if c <= tol
                break
            end
        end
        arg += log(fn_val)
        fn_val = zero
        if arg >= elimit
            fn_val = exp(arg)
        end
    else
        # Use a continued fraction expansion
        arg = p * log(x) - x - loggamma(p)
        a = one - p
        b = a + x + one
        c = zero
        pn1 = one
        pn2 = x
        pn3 = x + one
        pn4 = x * b
        fn_val = pn3 / pn4
        while true
            a += one
            b += two
            c += one
            an = a * c
            pn5 = b * pn3 - an * pn1
            pn6 = b * pn4 - an * pn2
            if abs(pn6) > zero
                rn = pn5 / pn6
                if abs(fn_val - rn) <= min(tol, tol * rn)
                    break
                end
                fn_val = rn
            end
            pn1 = pn3
            pn2 = pn4
            pn3 = pn5
            pn4 = pn6
            if abs(pn5) >= oflo
                # Re-scale terms in continued fraction if terms are large
                pn1 /= oflo
                pn2 /= oflo
                pn3 /= oflo
                pn4 /= oflo
            end
        end
        arg += log(fn_val)
        fn_val = one
        if arg >= elimit
            fn_val = one - exp(arg)
        end
    end

    return fn_val
end

end  # module GammaDistribution
