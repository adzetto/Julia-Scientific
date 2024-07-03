module RangeDistribution

export rngpi

using LinearAlgebra

function alnorm(x::Float64, upper::Bool)::Float64
    # Constants
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

function risf(x::Float64, t::Float64, n::Int)::Float64
    half = 0.5
    exp_factor = exp(-half * x * x)
    alnorm_factor = alnorm(x, false) - alnorm(x - t, false)
    fn_val = 0.3989422804 * exp_factor * (alnorm_factor)^(n - 1)
    return fn_val
end

function rngpi(t::Float64, n::Int)::Float64
    zero = 0.0
    half = 0.5
    two = 2.0
    eight = 8.0

    g = [0.4947004675, 0.4722875115, 0.4328156012, 0.3777022042, 0.3089381222, 0.2290083888, 0.1408017754, 0.04750625492]
    h = [0.01357622971, 0.03112676197, 0.04757925584, 0.06231448563, 0.07479799441, 0.08457825969, 0.09130170752, 0.09472530523]

    fn_val = zero
    if t <= zero || n <= 1
        return fn_val
    end

    xl = half * t
    a = half * (eight + xl)
    b = eight - xl
    y = zero

    for i in 1:8
        c = b * g[i]
        y += h[i] * (risf(a + c, t, n) + risf(a - c, t, n))
    end

    fn_val = (two * (alnorm(xl, false) - half))^n + two * b * y * n
    return fn_val
end

end # module RangeDistribution

using .RangeDistribution

# Example usage
t = 2.0
n = 10
fn_val = RangeDistribution.rngpi(t, n)
println("The probability of the normal range given t = $t and n = $n is: $fn_val")
