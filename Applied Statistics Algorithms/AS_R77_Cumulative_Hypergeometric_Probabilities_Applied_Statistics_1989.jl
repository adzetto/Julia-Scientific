module Hypergeometric

export chyper

function chyper(point::Bool, kk::Int, ll::Int, mm::Int, nn::Int)
    zero = 0.0
    half = 0.5
    one = 1.0
    sxteen = 16.0
    scale = 1.0e35
    rootpi = 2.506628274631001
    hundrd = 100.0
    elimit = -88.0
    mvbig = 1000
    mbig = 600

    k = kk + 1
    l = ll + 1
    m = mm + 1
    n = nn + 1
    dir = true

    # Check arguments are within permitted limits
    if n < 1 || m < n || k < 1 || k > m
        return zero, 1
    end

    if l < 1 || k - l > m - n
        return zero, 2
    end

    if !point
        return one, 0
    end

    if l > n || l > k
        return one, 2
    end

    if k == 1 || k == m || n == 1 || n == m
        return one, 0
    end

    if !point && ll == min(kk, nn)
        return one, 0
    end

    p = float(nn) / float(mm - nn)
    if float(min(kk, mm - kk)) > sxteen * max(p, one / p) && mm > mvbig && elimit > -hundrd
        # Use a normal approximation
        mean = float(kk) * float(nn) / float(mm)
        sig = sqrt(mean * (float(mm - nn) / float(mm)) * (float(mm - kk) / (float(mm - 1))))
        if point
            arg = -half * ((float(ll) - mean) / sig)^2
            fn_val = zero
            if arg >= elimit
                fn_val = exp(arg) / (sig * rootpi)
            end
        else
            fn_val = alnorm((float(ll) + half - mean) / sig, false)
        end
        return fn_val, 0
    else
        # Calculate exact hypergeometric probabilities
        if min(k - 1, m - k) > min(n - 1, m - n)
            i = k
            k = n
            n = i
        end

        if m - k < k - 1
            dir = !dir
            l = n - l + 1
            k = m - k + 1
        end

        if mm > mbig
            # Take logarithms of factorials
            p = alnfac(nn) - alnfac(mm) + alnfac(mm - kk) + alnfac(kk) +
                alnfac(mm - nn) - alnfac(ll) - alnfac(nn - ll) - alnfac(kk - ll) -
                alnfac(mm - nn - kk + ll)
            fn_val = zero
            if p >= elimit
                fn_val = exp(p)
            end
        else
            # Use Freeman/Lund algorithm
            fn_val = one
            for i in 1:(l - 1)
                fn_val *= float(k - i) * float(n - i) / (float(l - i) * float(m - i))
            end

            if l != k
                j = m - n + l
                for i in l:(k - 1)
                    fn_val *= float(j - i) / float(m - i)
                end
            end
        end

        if point
            return fn_val, 0
        end

        if fn_val == zero
            # We must recompute the point probability since it has underflowed
            if mm <= mbig
                p = alnfac(nn) - alnfac(mm) + alnfac(kk) + alnfac(mm - nn) -
                    alnfac(ll) - alnfac(nn - ll) - alnfac(kk - ll) - alnfac(mm - nn - kk + ll) +
                    alnfac(mm - kk)
            end
            p += log(scale)
            if p < elimit
                if ll > float(nn * kk + nn + kk + 1) / (mm + 2)
                    return one, 3
                end
                return zero, 3
            else
                p = exp(p)
            end
        else
            # Scale up at this point
            p = fn_val * scale
        end

        pt = zero
        nl = n - l
        kl = k - l
        mnkl = m - n - kl + 1

        if l <= kl
            for i in 1:(l - 1)
                p *= float(l - i) * float(mnkl - i) / (float(nl + i) * float(kl + i))
                pt += p
            end
            if p == zero
                return one, 3
            end
        else
            dir = !dir
            for j in 0:(kl - 1)
                p *= float(nl - j) * float(kl - j) / (float(l + j) * float(mnkl + j))
                pt += p
            end
            if p == zero
                return one, 3
            end
        end

        if dir
            fn_val += pt / scale
        else
            fn_val = one - pt / scale
        end

        return fn_val, 0
    end
end

function alnfac(n::Int)
    # This function calculates the logarithm of factorial n (ln(n!))
    if n <= 1
        return 0.0
    else
        return sum(log(i) for i in 2:n)
    end
end

function alnorm(x::Float64, upper::Bool)
    # This function evaluates the tail area of the standard normal curve from x to infinity
    # if upper is true or from minus infinity to x if upper is false.
    ltone = 7.0
    utzero = 18.66
    con = 1.28
    a1 = 0.398942280444
    a2 = 0.39990348504
    a3 = 5.75885480458
    a4 = 29.8213557807
    a5 = 2.62433121679
    b1 = 0.398942280385
    b2 = 3.8052e-8
    b3 = 1.00000615302
    b4 = 3.98064794e-4
    b5 = 1.98615381364
    b6 = 5.29330324926
    b7 = 0.151679116635
    b8 = 1.00000615302

    z = abs(x)
    if z <= ltone || (upper && z <= utzero)
        y = 0.5 * z * z
        if z > con
            alnorm = b1 * exp(-y) / (z + b2 + b3 / (z + b4 + b5 / (z + b6 + b7 / (z + b8))))
        else
            alnorm = 0.5 - z * (a1 - a2 * y / (y + a3 - a4 / (y + a5)))
        end
    else
        alnorm = 0.0
    end

    if !upper && x > 0.0
        alnorm = 1.0 - alnorm
    end

    if upper && x < 0.0
        alnorm = 1.0 - alnorm
    end

    return alnorm
end

end # module Hypergeometric

using .Hypergeometric

# Example usage
point = true
kk = 5
ll = 3
mm = 10
nn = 4
fn_val, ifault = Hypergeometric.chyper(point, kk, ll, mm, nn)

println("Function value: ", fn_val)
println("Ifault: ", ifault)
