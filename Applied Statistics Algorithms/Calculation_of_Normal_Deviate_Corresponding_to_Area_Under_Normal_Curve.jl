module NormalDeviate

using Printf, Random

# Subroutine to calculate the normal deviate for single precision
function ppnd7(p::Float32)::Tuple{Float32, Int}
    zero = 0.0f0
    one = 1.0f0
    half = 0.5f0
    split1 = 0.425f0
    split2 = 5.0f0
    const1 = 0.180625f0
    const2 = 1.6f0

    # Coefficients for P close to 0.5
    a = [3.3871327179f0, 5.0434271938f1, 1.5929113202f2, 5.9109374720f1]
    b = [1.7895169469f1, 7.8757757664f1, 6.7187563600f1]
    # Coefficients for P not close to 0, 0.5 or 1.
    c = [1.4234372777f0, 2.7568153900f0, 1.3067284816f0, 1.7023821103f-1]
    d = [7.3700164250f-1, 1.2021132975f-1]
    # Coefficients for P near 0 or 1.
    e = [6.6579051150f0, 3.0812263860f0, 4.2868294337f-1, 1.7337203997f-2]
    f = [2.4197894225f-1, 1.2258202635f-2]

    ifault = 0
    normal_dev = zero
    q = p - half

    if abs(q) <= split1
        r = const1 - q * q
        normal_dev = q * (((a[4] * r + a[3]) * r + a[2]) * r + a[1]) * r + a[1] / (((b[3] * r + b[2]) * r + b[1]) * r + one)
        return (normal_dev, ifault)
    else
        if q < zero
            r = p
        else
            r = one - p
        end

        if r <= zero
            ifault = 1
            return (zero, ifault)
        end

        r = sqrt(-log(r))

        if r <= split2
            r -= const2
            normal_dev = (((c[4] * r + c[3]) * r + c[2]) * r + c[1]) / ((d[2] * r + d[1]) * r + one)
        else
            r -= split2
            normal_dev = (((e[4] * r + e[3]) * r + e[2]) * r + e[1]) / ((f[2] * r + f[1]) * r + one)
        end

        if q < zero
            normal_dev = -normal_dev
        end

        return (normal_dev, ifault)
    end
end

# Subroutine to calculate the normal deviate for double precision
function ppnd16(p::Float64)::Tuple{Float64, Int}
    zero = 0.0
    one = 1.0
    half = 0.5
    split1 = 0.425
    split2 = 5.0
    const1 = 0.180625
    const2 = 1.6

    # Coefficients for P close to 0.5
    a = [3.3871328727963666080, 1.3314166789178437745e2, 1.9715909503065514427e3, 1.3731693765509461125e4,
         4.5921953931549871457e4, 6.7265770927008700853e4, 3.3430575583588128105e4, 2.5090809287301226727e3]
    b = [4.2313330701600911252e1, 6.8718700749205790830e2, 5.3941960214247511077e3, 2.1213794301586595867e4,
         3.9307895800092710610e4, 2.8729085735721942674e4, 5.2264952788528545610e3]
    # Coefficients for P not close to 0, 0.5 or 1.
    c = [1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550, 3.64784832476320460504,
         1.27045825245236838258, 2.41780725177450611770e-1, 2.27238449892691845833e-2, 7.74545014278341407640e-4]
    d = [2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1, 1.48103976427480074590e-1,
         1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9]
    # Coefficients for P near 0 or 1.
    e = [6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580, 2.96560571828504891230e-1,
         2.65321895265761230930e-2, 1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7]
    f = [5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2, 7.86869131145613259100e-4,
         1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15]

    ifault = 0
    normal_dev = zero
    q = p - half

    if abs(q) <= split1
        r = const1 - q * q
        normal_dev = q * (((((((a[8]*r + a[7])*r + a[6])*r + a[5])*r + a[4])*r + a[3])*r + a[2])*r + a[1]) * r + a[1] / (((((((b[7]*r + b[6])*r + b[5])*r + b[4])*r + b[3])*r + b[2])*r + b[1]) * r + one)
        return (normal_dev, ifault)
    else
        if q < zero
            r = p
        else
            r = one - p
        end

        if r <= zero
            ifault = 1
            return (zero, ifault)
        end

        r = sqrt(-log(r))

        if r <= split2
            r -= const2
            normal_dev = (((((((c[8]*r + c[7])*r + c[6])*r + c[5])*r + c[4])*r + c[3])*r + c[2])*r + c[1]) / (((((((d[8]*r + d[7])*r + d[6])*r + d[5])*r + d[4])*r + d[3])*r + d[2])*r + d[1]) * r + one)
        else
            r -= split2
            normal_dev = (((((((e[8]*r + e[7])*r + e[6])*r + e[5])*r + e[4])*r + e[3])*r + e[2])*r + e[1]) / (((((((f[8]*r + f[7])*r + f[6])*r + f[5])*r + f[4])*r + f[3])*r + f[2])*r + f[1]) * r + one)
        end

        if q < zero
            normal_dev = -normal_dev
        end

        return (normal_dev, ifault)
    end
end

end  # module NormalDeviate

# Example program to use the ppnd7 and ppnd16 functions
using .NormalDeviate

function main()
    while true
        println("Enter area under normal curve: ")
        p_dble = parse(Float64, readline())
        p_single = Float32(p_dble)
        
        z_single, ifault = NormalDeviate.ppnd7(p_single)
        if ifault != 0
            println("IFAULT = ", ifault)
            continue
        end
        
        z_dble, ifault = NormalDeviate.ppnd16(p_dble)
        if ifault != 0
            println("IFAULT = ", ifault)
            continue
        end
        
        @printf "Low prec. result = %11.6f  High prec. result = %16.11f\n" z_single z_dble
    end
end

main()
