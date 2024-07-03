module NormalScores

using Printf

# Function to calculate normal scores (Algorithm AS 177)
function nscor1(s::Vector{Float64}, n::Int, n2::Int, work::Array{Float64, 2})
    ifault = 3
    if n2 != div(n, 2)
        return s, ifault
    end

    ifault = 1
    if n <= 1
        return s, ifault
    end

    ifault = 0
    if n > 2000
        ifault = 2
    end

    an = n
    c = log(an)
    h = 0.025

    # Accumulate ordinates for calculation of integral for rankits
    for i in 1:n2
        i1 = i - 1
        ni = n - i
        ai1 = i1
        ani = ni
        scor = 0.0
        for j in 1:size(work, 2)
            scor += exp(work[2, j] + ai1 * work[3, j] + ani * work[4, j] + c) * work[1, j]
        end
        s[i] = scor * h
        c += log(ani / i)
    end

    return s, ifault
end

# Initialize the work array (Algorithm AS 177.1)
function init(work::Array{Float64, 2})
    xstart = -9.0
    h = 0.025
    pi2 = -0.918938533
    half = 0.5

    xx = xstart

    # Set up arrays for calculation of integral
    for i in 1:size(work, 2)
        work[1, i] = xx
        work[2, i] = pi2 - xx * xx * half
        work[3, i] = log(alnorm(xx, true))
        work[4, i] = log(alnorm(xx, false))
        xx = xstart + i * h
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

# Example usage
function main()
    n = 10
    n2 = div(n, 2)
    s = zeros(Float64, n2)
    work = zeros(Float64, 4, 721)
    
    # Initialize work array
    init(work)
    
    # Calculate normal scores
    s, ifault = nscor1(s, n, n2, work)
    
    # Print results
    println("Normal Scores: ", s)
    println("Ifault: ", ifault)
end

end

# Run the example
NormalScores.main()
