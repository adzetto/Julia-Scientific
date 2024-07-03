module Alnorm

export alnorm

function alnorm(x::Float64, upper::Bool)::Float64
    # Local variables
    zero = 0.0
    one = 1.0
    half = 0.5
    con = 1.28
    z = x
    y = 0.0
    up = upper
    fn_val = 0.0

    # Machine dependent constants
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

end # module Alnorm

using .Alnorm

x = 1.5
upper = true
fn_val = Alnorm.alnorm(x, upper)
println("The tail area of the standardised normal curve from x to infinity is: $fn_val")