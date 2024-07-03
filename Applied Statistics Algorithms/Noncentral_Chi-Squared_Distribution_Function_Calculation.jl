module Chi2NC

# Function to compute the non-central chi-square distribution function
function chi2nc(x::Float64, f::Float64, theta::Float64)::Float64
    # Local variables
    itrmax = 50
    errmax = 1.0e-12
    zero = 0.0
    one = 1.0
    two = 2.0

    fn_val = x

    if f <= zero || theta < zero
        println("Error in CHI2NC, IFAULT = 2")
        return fn_val
    end

    if x < zero
        println("Error in CHI2NC, IFAULT = 3")
        return fn_val
    end

    if x == zero
        return fn_val
    end

    lam = theta / two

    # Evaluate the first term
    n = one
    u = exp(-lam)
    v = u
    x2 = x / two
    f2 = f / two
    t = x2^f2 * exp(-x2) / exp(lngamma(f2 + one))
    term = v * t
    fn_val = term

    # Check if (f+2n) is greater than x
    flag = false
    while true
        if f + two * n - x > zero
            flag = true
            bound = t * x / (f + two * n - x)
            if bound <= errmax || n > itrmax
                if bound > errmax
                    println("Error in CHI2NC, IFAULT = 1")
                end
                return fn_val
            end
        end

        # Evaluate the next term of the expansion and then the partial sum
        u *= lam / n
        v += u
        t *= x / (f + two * n)
        term = v * t
        fn_val += term
        n += one

        if flag
            bound = t * x / (f + two * n - x)
            if bound <= errmax || n > itrmax
                if bound > errmax
                    println("Error in CHI2NC, IFAULT = 1")
                end
                return fn_val
            end
        end
    end
end

# Function to compute the logarithm of the gamma function using Lanczos approximation
function lngamma(z::Float64)::Float64
    # Coefficients for Lanczos approximation
    a = [0.9999999999995183, 676.5203681218835, -1259.139216722289, 771.3234287757674,
         -176.6150291498386, 12.50734324009056, -0.1385710331296526, 0.9934937113930748e-05,
         0.1659470187408462e-06]
    lnsqrt2pi = 0.9189385332046727
    half = 0.5
    sixpt5 = 6.5
    seven = 7.0

    if z <= 0.0
        println("Error: zero or -ve argument for lngamma")
        return 0.0
    end

    lanczos = 0.0
    tmp = z + seven
    for j in 9:-1:2
        lanczos += a[j] / tmp
        tmp -= 1.0
    end
    lanczos += a[1]
    lanczos = log(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half) * log(z + sixpt5)
    return lanczos
end

end  # module Chi2NC

# Example usage of the chi2nc function
using .Chi2NC

function main()
    x = 5.0
    f = 2.0
    theta = 3.0
    result = Chi2NC.chi2nc(x, f, theta)
    println("The value of the non-central chi-square distribution function is $result")
end

main()
