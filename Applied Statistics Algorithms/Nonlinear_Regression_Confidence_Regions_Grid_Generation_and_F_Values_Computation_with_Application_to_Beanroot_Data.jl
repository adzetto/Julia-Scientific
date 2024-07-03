module NonlinConfRegions

using LinearAlgebra

function halprn(param::Vector{Float64}, np::Int, ncols::Int, nobs::Int, np1::Int, np2::Int, p1lo::Float64, p1hi::Float64, p2lo::Float64, p2hi::Float64, conf::Vector{Float64}, nconf::Int, f::Vector{Float64}, grid::Matrix{Float64}, dim1::Int, ngrid1::Int, ngrid2::Int, deriv::Function, lindep::Vector{Bool}, ifault::Vector{Int})
    zero = 0.0
    one = 1.0
    hundrd = 100.0
    half = 0.5
    two = 2.0
    wt = 1.0

    ifault[1] = 0
    if ncols < 2
        ifault[1] = 1
    end
    if nobs <= ncols
        ifault[1] += 2
    end
    if np1 < 1 || np1 > np
        ifault[1] += 4
    end
    if np2 < 1 || np2 > np
        ifault[1] += 8
    end
    if np1 == np2
        ifault[1] += 16
    end
    if ngrid1 < 2 || ngrid2 < 2 || ngrid1 > dim1
        ifault[1] += 32
    end
    if ifault[1] > 0
        return
    end

    step1 = (p1hi - p1lo) / (ngrid1 - 1)
    step2 = (p2hi - p2lo) / (ngrid2 - 1)
    p1save = param[np1]
    p2save = param[np2]

    param[np1] = p1lo
    grad = zeros(Float64, ncols)
    for i1 in 1:ngrid1
        param[np2] = p2lo
        for i2 in 1:ngrid2
            # Initialize orthogonal reduction.
            include!()

            for iobs in 1:nobs
                deriv(param, np, ncols, iobs, grad, resid, wt)

                if np1 < ncols - 1
                    if np2 != ncols - 1
                        q = grad[ncols]
                        grad[ncols] = grad[np1]
                        grad[np1] = q
                    else
                        q = grad[ncols]
                        grad[ncols] = grad[np1]
                        grad[np1] = q
                    end
                end
                if np2 < ncols - 1
                    q = grad[ncols]
                    grad[ncols] = grad[np2]
                    grad[np2] = q
                end

                include!(wt, grad, resid)
            end

            tolset(1.0e-8)
            sing(lindep, ifault)

            ndf = nobs - ncols - ifault[1]

            ssq1 = if ncols > 2
                rss[ncols - 2]
            else
                rss[2]
            end

            grid[i1, i2] = half * (ssq1 - sserr) / (sserr / ndf)
            param[np2] += step2
        end
        param[np1] += step1
    end

    param[np1] = p1save
    param[np2] = p2save

    for i in 1:nconf
        q = one - conf[i] / hundrd
        if q <= zero || q > one
            ifault[1] = -i
            return
        end
        f[i] = half * ndf * (q^(-two / ndf) - one)
    end

    return
end

end

module MyData

using Random

const x = zeros(Float64, 100)
const y = zeros(Float64, 100)

end

using .NonlinConfRegions
using .MyData

function logistic4(param::Vector{Float64}, np::Int, ncols::Int, iobs::Int, grad::Vector{Float64}, resid::Float64, wt::Float64)
    diff = MyData.x[iobs] - param[4]
    expntl = exp(-param[3] * diff)
    denom = 1.0 + expntl
    resid = MyData.y[iobs] - param[1] - (param[2] - param[1]) / denom
    grad[1] = -1.0 + 1.0 / denom
    grad[2] = -1.0 / denom
    grad[3] = -(param[2] - param[1]) * diff * expntl / denom^2
    grad[4] = (param[2] - param[1]) * param[3] * expntl / denom^2
end

function main()
    fname = "beanroot.dat"
    fstatus = 0
    nobs = 0
    open(fname, "r") do file
        for line in eachline(file)
            MyData.x[nobs + 1], MyData.y[nobs + 1] = split(line) .|> parse
            nobs += 1
        end
    end
    println("No. of cases read = ", nobs)

    param = [0.88, 21.3, 0.0, 0.0]
    np = 4
    ncols = 4
    np1 = 3
    np2 = 4
    p1lo = 0.4
    p1hi = 1.5
    p2lo = 5.8
    p2hi = 7.2
    conf = [10.0, 50.0, 90.0, 95.0, 98.0, 99.0]
    nconf = 6
    dim1 = 31
    ngrid1 = 31
    ngrid2 = 31
    grid = zeros(Float64, dim1, ngrid2)
    f = zeros(Float64, nconf)
    lindep = falses(ncols)
    ifault = zeros(Int, 1)

    halprn(param, np, ncols, nobs, np1, np2, p1lo, p1hi, p2lo, p2hi, conf, nconf, f, grid, dim1, ngrid1, ngrid2, logistic4, lindep, ifault)

    println("Percentage points of the F-distribution")
    println("   %    F-value")
    for i in 1:nconf
        println(conf[i], "    ", f[i])
    end

    open("confgrid.dat", "w") do file
        for row in 1:ngrid1
            println(file, join(grid[row, :], " "))
        end
    end
end

main()
