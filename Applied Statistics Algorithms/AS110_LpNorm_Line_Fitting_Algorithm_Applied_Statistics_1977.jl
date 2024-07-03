module LpNormFit

export lpest

function lpest(n::Int, p::Float64, x::Vector{Float64}, y::Vector{Float64}, maxit::Int)
    # Outputs
    a = 0.0
    b = 0.0
    sd = 0.0
    r = ones(Float64, n)
    rate = 0.0
    it = 0
    npo = 0
    ifault = 0

    # Constants
    eps = 1.0e-6
    wp = p - 2.0
    eps2 = 2.0 * eps

    # Check input validity
    if maxit < 2
        ifault = 2
        return a, b, sd, r, rate, it, npo, ifault
    end

    if n < 2
        ifault = 4
        return a, b, sd, r, rate, it, npo, ifault
    end

    for it in 1:maxit
        npo = 0
        sw = 0.0
        xmean = 0.0
        ymean = 0.0
        ssx = 0.0
        spxy = 0.0

        for i in 1:n
            absri = abs(r[i])
            if absri <= eps
                npo += 1
                continue
            end

            w = absri^wp
            sw += w
            div = w / sw
            xi = x[i] - xmean
            yi = y[i] - ymean
            xiw = xi * w
            dx = xi * xiw
            dxy = yi * xiw
            ssx += dx - dx * div
            spxy += dxy - dxy * div
            xmean += xi * div
            ymean += yi * div
        end

        if ssx < eps
            ifault = 3
            return a, b, sd, r, rate, it, npo, ifault
        end

        b = spxy / ssx
        a = ymean - b * xmean

        sd2 = 0.0
        isw = 0

        for i in 1:n
            res = y[i] - a - b * x[i]
            absri = abs(res)
            if abs(absri - abs(r[i])) > eps2
                isw = 1
            end
            sd2 += absri^p
            r[i] = res
        end

        rate = abs(sd2 - sd) / sd2

        if isw == 0
            return a, b, sd, r, rate, it, npo, ifault
        end

        if it == 1
            sd = sd2
            a2 = a
            b2 = b
            continue
        end

        if sd2 > sd
            ifault = 1
            a = a2
            b = b2
            r .= y .- a .- b .* x
            return a, b, sd, r, rate, it, npo, ifault
        end

        sd = sd2
        a2 = a
        b2 = b
    end

    ifault = 5
    return a, b, sd, r, rate, it, npo, ifault
end

end # module LpNormFit

# Usage example
using .LpNormFit

n = 100
p = 1.5
x = randn(n)
y = 2.0 .+ 3.0 .* x .+ 0.5 .* randn(n)
maxit = 100

a, b, sd, r, rate, it, npo, ifault = LpNormFit.lpest(n, p, x, y, maxit)

println("Estimates:")
println("a = ", a)
println("b = ", b)
println("sd = ", sd)
println("rate = ", rate)
println("iterations = ", it)
println("npo = ", npo)
println("ifault = ", ifault)
