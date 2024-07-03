# Load necessary packages
using LinearAlgebra

# Module for Errors in Variables
module ErrorsInVariables

export evms, bf, zed

# User defined functions bf and zed
function bf(b, f, theta, xi)
    dp = Float64

    # The parameters in theta are:
    # theta[1] = x0
    # theta[2] = y0
    # theta[3] = r
    # theta[4] = angle 1
    # theta[5] = angle 2

    # The (x,y,z) coordinates are xi[1], xi[2] & xi[3]

    two = 2.0

    # Calculate F
    sine1 = sin(theta[4])
    cos1  = cos(theta[4])
    sine2 = sin(theta[5])
    cos2  = cos(theta[5])
    x = xi[1] - theta[1]
    y = xi[2] - theta[2]
    z = xi[3]
    e1 = x*cos1 + (z*cos2 - y*sine2)*sine1
    e2 = y*cos2 + z*sine2
    f[1] = e1^2 + e2^2 - theta[3]^2

    # Calculate derivatives w.r.t. x, y, z
    b[1,1] = two*e1*cos1
    b[1,2] = two*e2*cos2 - two*e1*sine2*sine1
    b[1,3] = two*e1*cos2*sine1 + two*e2*sine2

    return
end

function zed(b, f, theta, xi, grad)
    dp = Float64

    two = 2.0

    sine1 = sin(theta[4])
    cos1  = cos(theta[4])
    sine2 = sin(theta[5])
    cos2  = cos(theta[5])
    x = xi[1] - theta[1]
    y = xi[2] - theta[2]
    z = xi[3]
    e1 = x*cos1 + (z*cos2 - y*sine2)*sine1
    e2 = y*cos2 + z*sine2

    # Calculate derivatives w.r.t. parameters
    grad[1,1] = -two*e1*cos1
    grad[1,2] = -two*e2*cos2 + two*e1*sine2*sine1
    grad[1,3] = -two*theta[3]
    grad[1,4] = two*e1*(-x*sine1 + (z*cos2 - y*sine2)*cos1)
    grad[1,5] = two*e2*(-y*sine2 + z*cos2) + two*e1*(-z*sine2 - y*cos2)*sine1

    return
end

# EIV Module containing the core algorithm
function evms(bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, theta, v)

    crit1 = 1.e-6
    one = 1.0
    zero = 0.0
    niter1 = 50

    ifault[1] = 0
    ifault[2] = 0

    for iter in 1:niter1
        inner(bs, DATA, est, g1, ifault, ndat, npar, nvar, resid, theta, v)

        g2 = copy(g1)
        q2 = copy(g1[:, 1])

        chol!(g2)
        backsub!(g2, q2, 1, npar, 1)

        wdiff = zero
        for ipar in 1:npar
            w1 = q2[ipar]
            w2 = theta[ipar]
            wdiff = max(wdiff, abs(w1 / w2))
            theta[ipar] = w2 - w1
        end

        if wdiff <= crit1
            break
        end
    end

    ifault[1] = 1

    for ipar1 in 1:npar
        for ipar2 in ipar1:npar
            g2[ipar1, ipar2] = g1[ipar1, ipar2]
            g3[ipar1, ipar2] = zero
        end
    end

    chol!(g2)
    for ipar1 in 1:npar
        g3[ipar1, ipar1] = one / g2[ipar1, ipar1]
    end
    backsub!(g2, g3, npar, npar, 3)

    return
end

function inner(bs, DATA, est, g1, ifault, ndat, npar, nvar, resid, theta, v)
    crit2 = 1.e-6
    half = 0.5
    zero = 0.0
    niter2 = 20

    phi = zero
    q1 = zeros(Float64, npar)

    for ipar1 in 1:npar
        for ipar2 in ipar1:npar
            g1[ipar1, ipar2] = zero
        end
    end

    for idat in 1:ndat
        xi = DATA[idat, 1:nvar]
        for iter in 1:niter2
            wdiff = zero
            f = zeros(Float64, 1)
            b = zeros(Float64, 1, nvar)
            bf(b, f, theta, xi)

            bv = b * v
            s = bv * transpose(b)

            for ivar1 in 1:nvar
                h = f[1] + dot(b[1, 1:nvar], DATA[idat, 1:nvar] - xi)
                xi_new = DATA[idat, ivar1] - dot(bv[:, ivar1], h) / s[1, 1]
                wdiff = max(wdiff, abs((xi_new - xi[ivar1]) / xi_new))
                xi[ivar1] = xi_new
            end

            if wdiff <= crit2
                break
            end
        end

        ifault[2] = 1

        phi += h^2 / s[1, 1]

        est[idat, 1:nvar] .= xi
        resid[idat, 1:nvar] .= DATA[idat, 1:nvar] - xi

        z = zeros(Float64, 1, npar)
        zed(b, f, theta, xi, z)

        for ipar1 in 1:npar
            q1[ipar1] += z[1, ipar1] * h / s[1, 1]
            for ipar2 in ipar1:npar
                g1[ipar1, ipar2] += z[1, ipar1] * z[1, ipar2] / s[1, 1]
            end
        end
    end

    phi *= half

    return
end

function chol!(a)
    n = size(a, 1)
    for j in 1:n
        for i in 1:j-1
            a[j,j] -= a[i,j]^2
        end
        a[j,j] = sqrt(a[j,j])
        for k in j+1:n
            for i in 1:j-1
                a[j,k] -= a[i,j]*a[i,k]
            end
            a[j,k] /= a[j,j]
        end
    end
end

function backsub!(a, b, m, n, job)
    for ic1 in m:-1:1
        if job == 2
            k1, k2, k3 = 1, n, 1
        else
            k1, k2, k3 = n, 1, -1
            if job == 1
                k1 = n
            else
                k1 = ic1
            end
        end
        for ir1 in k1:k2
            w1 = b[ir1, ic1]
            if job == 2
                k2 = ir1 - 1
            else
                k1 = n
                k2 = ir1 + 1
            end
            for ir2 in k1:k2:k3
                w1 -= a[ir1, ir2] * b[ir2, ic1]
            end
            b[ir1, ic1] = w1 / a[ir1, ir1]
            if job == 3
                b[ic1, ir1] = b[ir1, ic1]
            end
        end
    end
end

end # module ErrorsInVariables

# Main Program for the Cylinder Example
function main_cylinder()
    using Random
    using LinearAlgebra
    using ErrorsInVariables

    Random.seed!(1234)  # For reproducibility

    ndat = 50
    npar = 5
    nvar = 3
    ifault = zeros(Int, 2)
    bs = zeros(Float64, 1, nvar)
    DATA = zeros(Float64, ndat, nvar)
    est = zeros(Float64, ndat, nvar)
    g1 = zeros(Float64, npar, npar)
    g3 = zeros(Float64, npar, npar)
    theta = [1.5, 2.5, 3.5, 1.0, 1.0]
    resid = zeros(Float64, ndat, nvar)
    v = Matrix{Float64}(I, nvar, nvar)

    # True parameter values: 1.0, 2.0, 3.0, 0.6, 1.2
    sine1 = sin(0.6)
    cos1  = cos(0.6)
    sine2 = sin(1.2)
    cos2  = cos(1.2)

    for i in 1:ndat
        e1 = rand() - 0.5
        e2 = rand() - 0.5
        scale = sqrt(9.0 / (e1^2 + e2^2))
        e1 *= scale
        e2 *= scale
        z = 5.0 * rand()
        y = (e2 - z * sine2) / cos2
        x = (e1 - (z * cos2 - y * sine2) * sine1) / cos1
        DATA[i, :] = [x + 1.0, y + 2.0, z]
    end

    # Perturb the data slightly
    for i in 1:ndat
        DATA[i, :] .+= 0.25 .* (rand(3) .- 0.5)
    end

    # Call the EVMS function
    ErrorsInVariables.evms(bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, theta, v)

    println("Error status: ", ifault)
    println("Parameter estimates: ", theta)
    println("Information matrix:")
    println(g1)
    println("Inverse of Information matrix:")
    println(g3)
    println()
    println("   X      Y      Z           Estimates              Residuals")
    for i in 1:ndat
        println(DATA[i, :], "  ", est[i, :], "  ", resid[i, :])
    end
end

main_cylinder()


using LinearAlgebra
using Random

# Module for Errors in Variables
module ErrorsInVariables

export evms, evm, bf, zed

function bf_cylinder(b, f, theta, xi)
    two = 2.0
    sine1 = sin(theta[4])
    cos1  = cos(theta[4])
    sine2 = sin(theta[5])
    cos2  = cos(theta[5])
    x = xi[1] - theta[1]
    y = xi[2] - theta[2]
    z = xi[3]
    e1 = x*cos1 + (z*cos2 - y*sine2)*sine1
    e2 = y*cos2 + z*sine2
    f[1] = e1^2 + e2^2 - theta[3]^2

    b[1,1] = two*e1*cos1
    b[1,2] = two*e2*cos2 - two*e1*sine2*sine1
    b[1,3] = two*e1*cos2*sine1 + two*e2*sine2
end

function zed_cylinder(b, f, theta, xi, grad)
    two = 2.0
    sine1 = sin(theta[4])
    cos1  = cos(theta[4])
    sine2 = sin(theta[5])
    cos2  = cos(theta[5])
    x = xi[1] - theta[1]
    y = xi[2] - theta[2]
    z = xi[3]
    e1 = x*cos1 + (z*cos2 - y*sine2)*sine1
    e2 = y*cos2 + z*sine2

    grad[1,1] = -two*e1*cos1
    grad[1,2] = -two*e2*cos2 + two*e1*sine2*sine1
    grad[1,3] = -two*theta[3]
    grad[1,4] = two*e1*(-x*sine1 + (z*cos2 - y*sine2)*cos1)
    grad[1,5] = two*e2*(-y*sine2 + z*cos2) + two*e1*(-z*sine2 - y*cos2)*sine1
end

function bf_ex1(b, f, theta, xi)
    f[1]=theta[1]^2 + theta[2]*xi[1] + theta[3]*xi[2]^2 - xi[3]
    f[2]=theta[2]^2 + (theta[3]*xi[1])^2 + theta[4]^3*xi[3]^2 - xi[4]
    f[3]=theta[3]^2 + (theta[1]*xi[2])^2 - xi[5]
    b[1,1]=theta[2]
    b[1,2]=2.0*theta[3]*xi[2]
    b[1,3]=-1.0
    b[1,4]=0.0
    b[1,5]=0.0
    b[2,1]=2.0*theta[3]^2*xi[1]
    b[2,2]=0.0
    b[2,3]=2.0*theta[4]^3*xi[3]
    b[2,4]=-1.0
    b[2,5]=0.0
    b[3,1]=0.0
    b[3,2]=2.0*theta[1]^2*xi[2]
    b[3,3]=0.0
    b[3,4]=0.0
    b[3,5]=-1.0
end

function zed_ex1(b, f, theta, xi, z)
    z[1,1]=2.0*theta[1]
    z[1,2]=xi[1]
    z[1,3]=xi[2]^2
    z[1,4]=0.0
    z[2,1]=0.0
    z[2,2]=2.0*theta[2]
    z[2,3]=2.0*theta[3]*xi[1]^2
    z[2,4]=3.0*(theta[4]*xi[3])^2
    z[3,1]=2.0*theta[1]*xi[2]^2
    z[3,2]=0.0
    z[3,3]=2.0*theta[3]
    z[3,4]=0.0
end

function bf_ex2(b, f, theta, xi)
    w1=xi[1] - theta[1]
    w2=xi[2] - theta[2]
    f[1]=theta[3]*w1*w1 + 2.0*theta[4]*w1*w2 + theta[5]*w2*w2 - 1.0
    b[1,1]=2.0*(w1*theta[3] + w2*theta[4])
    b[1,2]=2.0*(w1*theta[4] + w2*theta[5])
end

function zed_ex2(b, f, theta, xi, z)
    w1=xi[1] - theta[1]
    w2=xi[2] - theta[2]
    z[1,1]=-2.0*(w1*theta[3] + w2*theta[4])
    z[1,2]=-2.0*(w1*theta[4] + w2*theta[5])
    z[1,3]=w1*w1
    z[1,4]=2.0*w1*w2
    z[1,5]=w2*w2
end

function evms(bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, theta, v)
    crit1 = 1.e-6
    one = 1.0
    zero = 0.0
    niter1 = 50

    ifault[1] = 0
    ifault[2] = 0

    for iter in 1:niter1
        inners(bs, DATA, est, g1, ifault, ndat, npar, nvar, resid, theta, v)

        g2 = copy(g1)
        q2 = copy(g1[:, 1])

        chol!(g2)
        backsub!(g2, q2, 1, npar, 1)

        wdiff = zero
        for ipar in 1:npar
            w1 = q2[ipar]
            w2 = theta[ipar]
            wdiff = max(wdiff, abs(w1 / w2))
            theta[ipar] = w2 - w1
        end

        if wdiff <= crit1
            break
        end
    end

    ifault[1] = 1

    for ipar1 in 1:npar
        for ipar2 in ipar1:npar
            g2[ipar1, ipar2] = g1[ipar1, ipar2]
            g3[ipar1, ipar2] = zero
        end
    end

    chol!(g2)
    for ipar1 in 1:npar
        g3[ipar1, ipar1] = one / g2[ipar1, ipar1]
    end
    backsub!(g2, g3, npar, npar, 3)

    return
end

function inners(bs, DATA, est, g1, ifault, ndat, npar, nvar, resid, theta, v)
    crit2 = 1.e-5
    half = 0.5
    zero = 0.0
    niter2 = 20

    phi = zero
    q1 = zeros(Float64, npar)

    for ipar1 in 1:npar
        for ipar2 in ipar1:npar
            g1[ipar1, ipar2] = zero
        end
    end

    for idat in 1:ndat
        xi = DATA[idat, 1:nvar]
        for iter in 1:niter2
            wdiff = zero
            f = zeros(Float64, 1)
            b = zeros(Float64, 1, nvar)
            bf(bs, f, theta, xi)

            bv = b * v
            s = bv * transpose(b)

            for ivar1 in 1:nvar
                h = f[1] + dot(b[1, 1:nvar], DATA[idat, 1:nvar] - xi)
                xi_new = DATA[idat, ivar1] - dot(bv[:, ivar1], h) / s[1, 1]
                wdiff = max(wdiff, abs((xi_new - xi[ivar1]) / xi_new))
                xi[ivar1] = xi_new
            end

            if wdiff <= crit2
                break
            end
        end

        ifault[2] = 1

        phi += f[1]^2 / s[1, 1]

        est[idat, 1:nvar] = xi
        resid[idat, 1:nvar] = DATA[idat, 1:nvar] - xi

        z = zeros(Float64, 1, npar)
        zed(bs, f, theta, xi, z)
        backsub!(s, z, 1, 1, 2)

        for ipar1 in 1:npar
            q1[ipar1] += z[1, ipar1] * f[1] / s[1, 1]
            for ipar2 in ipar1:npar
                g1[ipar1, ipar2] += z[1, ipar1] * z[1, ipar2] / s[1, 1]
            end
        end
    end

    phi *= half

    return
end

function chol!(a)
    n = size(a, 1)
    for j in 1:n
        a[j,j] = sqrt(a[j,j] - sum(a[1:j-1, j].^2))
        for i in j+1:n
            a[j,i] = (a[j,i] - sum(a[1:j-1, j] .* a[1:j-1, i])) / a[j,j]
        end
    end
    return a
end

function backsub!(a, b, m, n, job)
    for ic1 in m:-1:1
        if job == 2
            k1 = 1
            k2 = n
            k3 = 1
        else
            k2 = 1
            k3 = -1
            k1 = job == 1 ? n : ic1
        end
        for ir1 in k1:k2:k3
            w1 = b[ir1, ic1]
            if job == 2
                k2 = ir1 - 1
            else
                k1 = n
                k2 = ir1 + 1
            end
            for ir2 in k1:k2:k3
                w1 -= a[ir1, ir2] * b[ir2, ic1]
            end
            b[ir1, ic1] = w1 / a[ir1, ir1]
            if job == 3
                b[ic1, ir1] = b[ir1, ic1]
            end
        end
    end
    return b
end

end # End of module ErrorsInVariables

using .ErrorsInVariables

# Main function for Cylinder example
function main_cylinder()
    ndat = 50
    npar = 5
    nvar = 3
    ifault = zeros(Int, 2)
    bs = zeros(Float64, 1, nvar)
    DATA = zeros(Float64, ndat, nvar)
    est = zeros(Float64, ndat, nvar)
    g1 = zeros(Float64, npar, npar)
    g3 = zeros(Float64, npar, npar)
    theta = [1.5, 2.5, 3.5, 1.0, 1.0]
    resid = zeros(Float64, ndat, nvar)
    v = Matrix{Float64}(I, nvar, nvar)

    # True parameter values: 1.0, 2.0, 3.0, 0.6, 1.2
    sine1 = sin(0.6)
    cos1  = cos(0.6)
    sine2 = sin(1.2)
    cos2  = cos(1.2)

    for i in 1:ndat
        e1 = rand() - 0.5
        e2 = rand() - 0.5
        scale = sqrt(9.0 / (e1^2 + e2^2))
        e1 *= scale
        e2 *= scale
        z = 5.0 * rand()
        y = (e2 - z * sine2) / cos2
        x = (e1 - (z * cos2 - y * sine2) * sine1) / cos1
        DATA[i, :] = [x + 1.0, y + 2.0, z]
    end

    # Perturb the data slightly
    for i in 1:ndat
        DATA[i, :] .+= 0.25 .* (rand(3) .- 0.5)
    end

    # Call the EVMS function
    ErrorsInVariables.evms(bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, theta, v)

    println("Error status: ", ifault)
    println("Parameter estimates: ", theta)
    println("Information matrix:")
    println(g1)
    println("Inverse of Information matrix:")
    println(g3)
    println()
    println("   X      Y      Z           Estimates              Residuals")
    for i in 1:ndat
        println(DATA[i, :], "  ", est[i, :], "  ", resid[i, :])
    end
end

# Main function for EX1
function main_ex1()
    ndat = 6
    neq = 3
    npar = 4
    nvar = 5
    ifault = zeros(Int, 2)
    b = zeros(Float64, neq, nvar)
    DATA = [
        -0.78  1.39  3.99 11.70  3.44;
        -1.27  2.00  5.15 36.58  4.59;
        2.62  0.38  5.77 16.49  0.33;
        2.05  2.92  6.78 75.68  2.74;
        1.05  2.46 12.85 149.42  7.59;
        2.96  1.47  3.90 52.11  6.82
    ]
    est = zeros(Float64, ndat, nvar)
    g1 = zeros(Float64, npar, npar)
    g3 = zeros(Float64, npar, npar)
    theta = [1.0, 1.0, 1.0, 1.0]
    resid = zeros(Float64, ndat, nvar)
    v = [
        1.00  -0.20   0.10   0.30   0.00;
        -0.20   1.00   0.20  -0.10   0.00;
        0.10   0.20   2.00   0.40  -0.20;
        0.30  -0.10   0.40   5.00   0.40;
        0.00   0.00  -0.20   0.40   3.00
    ]

    ErrorsInVariables.evm(b, DATA, est, g1, g3, ifault, ndat, neq, npar, nvar, 0.0, zeros(Float64, npar), resid, theta, v)

    println("Error status: ", ifault)
    println("Parameter estimates: ", theta)
    println("Information matrix:")
    println(g1)
    println("Inverse of Information matrix:")
    println(g3)
    println()
    println("   DATA           ESTIMATES              RESIDUALS")
    for i in 1:ndat
        println(DATA[i, :], "  ", est[i, :], "  ", resid[i, :])
    end
end

# Main function for EX2
function main_ex2()
    ndat = 32
    npar = 5
    nvar = 2
    ifault = zeros(Int, 2)
    bs = zeros(Float64, 1, nvar)
    DATA = [
        1420.0 470.0;
        1188.0 476.0;
        1015.0 539.0;
        830.0 655.0;
        660.0 825.0;
        543.0 1041.0;
        468.0 1268.0;
        459.0 1528.0;
        510.0 1804.0;
        628.0 2063.0;
        795.0 2221.0;
        1040.0 2345.0;
        1303.0 2417.0;
        1452.0 2431.0;
        1738.0 2411.0;
        1901.0 2318.0;
        2091.0 2222.0;
        2230.0 2117.0;
        2311.0 1984.0;
        2385.0 1837.0;
        2451.0 1630.0;
        2467.0 1462.0;
        2472.0 1263.0;
        2446.0 1112.0;
        2387.0 963.0;
        2255.0 796.0;
        2109.0 696.0;
        1929.0 597.0;
        1846.0 556.0;
        1620.0 486.0;
        1418.0 437.0;
        1259.0 453.0
    ]
    est = zeros(Float64, ndat, nvar)
    g1 = zeros(Float64, npar, npar)
    g3 = zeros(Float64, npar, npar)
    theta = [1.0, 1.0, 1.0, 1.0, 1.0]
    resid = zeros(Float64, ndat, nvar)
    v = Matrix{Float64}(I, nvar, nvar)

    for i in 1:ndat
        DATA[i, :] .= 0.001 .* DATA[i, :]
    end

    ErrorsInVariables.evms(bs, DATA, est, g1, g3, ifault, ndat, npar, nvar, resid, theta, v)

    for i in 1:ndat
        DATA[i, :] .= 1000.0 .* DATA[i, :]
        est[i, :] .= 1000.0 .* est[i, :]
        resid[i, :] .= 1000.0 .* resid[i, :]
    end

    println("Error status: ", ifault)
    println("Parameter estimates: ", theta)
    println("Information matrix:")
    println(g1)
    println("Inverse of Information matrix:")
    println(g3)
    println()
    println("   DATA           ESTIMATES              RESIDUALS")
    for i in 1:ndat
        println(DATA[i, :], "  ", est[i, :], "  ", resid[i, :])
    end
end

# main_cylinder()
# main_ex1()
# main_ex2()
