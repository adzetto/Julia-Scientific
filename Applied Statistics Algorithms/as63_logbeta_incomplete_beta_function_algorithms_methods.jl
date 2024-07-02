# Translated from FORTRAN by adzetto

module LogBeta
# Log of the beta function
# Includes log of the gamma function

using Base.MathConstants: π

export betaln, betain

dp = Float64

function betaln(a0::dp, b0::dp)::dp
    #-----------------------------------------------------------------------
    #     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
    #-----------------------------------------------------------------------
    #     E = 0.5*LN(2*PI)
    #--------------------------

    # Local variables
    e = 0.5 * log(2π)
    a = min(a0, b0)
    b = max(a0, b0)
    
    if a < 8.0
        if a < 1.0
            #-----------------------------------------------------------------------
            #                   PROCEDURE WHEN A < 1
            #-----------------------------------------------------------------------
            if b < 8.0
                return gamln(a) + (gamln(b) - gamln(a + b))
            end
            return gamln(a) + algdiv(a, b)
        end
        #-----------------------------------------------------------------------
        #                PROCEDURE WHEN 1 <= A < 8
        #-----------------------------------------------------------------------
        if a <= 2.0
            if b <= 2.0
                return gamln(a) + gamln(b) - gsumln(a, b)
            end
            w = 0.0
            if b < 8.0
                n = a - 1.0
                w = 1.0
                for i in 1:Int(n)
                    a -= 1.0
                    h = a / b
                    w *= h / (1.0 + h)
                end
                w = log(w)
                if b >= 8.0
                    return w + gamln(a) + algdiv(a, b)
                end
                z = 1.0
                n = b - 1.0
                for i in 1:Int(n)
                    b -= 1.0
                    z *= b / (a + b)
                end
                return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)))
            end
            return gamln(a) + algdiv(a, b)
        end
        #                REDUCTION OF A WHEN B <= 1000

        if b > 1000.0
            n = a - 1.0
            w = 1.0
            for i in 1:Int(n)
                a -= 1.0
                w *= a / (1.0 + a / b)
            end
            return (log(w) - n * log(b)) + (gamln(a) + algdiv(a, b))
        end
    end
    #-----------------------------------------------------------------------
    #                   PROCEDURE WHEN A >= 8
    #-----------------------------------------------------------------------
    w = bcorr(a, b)
    h = a / b
    c = h / (1.0 + h)
    u = -(a - 0.5) * log(c)
    v = b * alnrel(h)
    if u > v
        return ((-0.5 * log(b) + e) + w - v) - u
    end
    return ((-0.5 * log(b) + e) + w - u) - v
end

function gamln(a::dp)::dp
    #-----------------------------------------------------------------------
    #      EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
    #-----------------------------------------------------------------------
    #            WRITTEN BY ALFRED H. MORRIS
    #            NAVAL SURFACE WARFARE CENTER
    #               DAHLGREN, VIRGINIA
    #--------------------------
    #              D = 0.5*(LN(2*PI) - 1)
    #--------------------------
    
    d = 0.5 * (log(2π) - 1)
    c = [0.833333333333333e-01, -0.277777777760991e-02, 0.793650666825390e-03,
               -0.595202931351870e-03, 0.837308034031215e-03, -0.165322962780713e-02]
    
    if a <= 0.8
        return gamln1(a) - log(a)
    end
    if a <= 2.25
        t = (a - 0.5) - 0.5
        return gamln1(t)
    end
    
    if a < 10.0
        n = a - 1.25
        t = a
        w = 1.0
        for i in 1:Int(n)
            t -= 1.0
            w *= t
        end
        return gamln1(t - 1.0) + log(w)
    end
    
    t = (1.0 / a) ^ 2
    w = ((((c[6] * t + c[5]) * t + c[4]) * t + c[3]) * t + c[2]) * t + c[1]
    return (d + w) + (a - 0.5) * (log(a) - 1.0)
end

function algdiv(a::dp, b::dp)::dp
    #-----------------------------------------------------------------------
    #     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8
    #
    #                         --------
    #
    #     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
    #     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
    #-----------------------------------------------------------------------
    
    c = [0.833333333333333e-01, -0.277777777760991e-02, 0.793650666825390e-03,
               -0.595202931351870e-03, 0.837308034031215e-03, -0.165322962780713e-02]
    
    if a > b
        h = b / a
        c = 1.0 / (1.0 + h)
        x = h / (1.0 + h)
        d = a + (b - 0.5)
    else
        h = a / b
        c = h / (1.0 + h)
        x = 1.0 / (1.0 + h)
        d = b + (a - 0.5)
    end
    
    x2 = x * x
    s3 = 1.0 + (x + x2)
    s5 = 1.0 + (x + x2 * s3)
    s7 = 1.0 + (x + x2 * s5)
    s9 = 1.0 + (x + x2 * s7)
    s11 = 1.0 + (x + x2 * s9)
    
    t = (1.0 / b) ^ 2
    w = ((((c[6] * s11 * t + c[5] * s9) * t + c[4] * s7) * t + c[3] * s5) * t + c[2] * s3) * t + c[1]
    w = w * (c / b)
    
    u = d * alnrel(a / b)
    v = a * (log(b) - 1.0)
    if u > v
        return (w - v) - u
    end
    return (w - u) - v
end

function gsumln(a::dp, b::dp)::dp
    #-----------------------------------------------------------------------
    #          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
    #          FOR 1 <= A <= 2  AND  1 <= B <= 2
    #-----------------------------------------------------------------------
    
    x = a + b - 2.0
    if x <= 0.25
        return gamln1(1.0 + x)
    end
    if x <= 1.25
        return gamln1(x) + alnrel(x)
    end
    return gamln1(x - 1.0) + log(x * (1.0 + x))
end

function bcorr(a0::dp, b0::dp)::dp
    #-----------------------------------------------------------------------
    #     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
    #     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
    #     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8.
    #-----------------------------------------------------------------------
    
    c = [0.833333333333333e-01, -0.277777777760991e-02, 0.793650666825390e-03,
               -0.595202931351870e-03, 0.837308034031215e-03, -0.165322962780713e-02]
    
    a = min(a0, b0)
    b = max(a0, b0)
    
    h = a / b
    cc = h / (1.0 + h)
    x = 1.0 / (1.0 + h)
    x2 = x * x
    
    s3 = 1.0 + (x + x2)
    s5 = 1.0 + (x + x2 * s3)
    s7 = 1.0 + (x + x2 * s5)
    s9 = 1.0 + (x + x2 * s7)
    s11 = 1.0 + (x + x2 * s9)
    
    t = (1.0 / b) ^ 2
    w = ((((c[6] * s11 * t + c[5] * s9) * t + c[4] * s7) * t + c[3] * s5) * t + c[2] * s3) * t + c[1]
    w = w * (cc / b)
    
    t = (1.0 / a) ^ 2
    return (((((c[6] * t + c[5]) * t + c[4]) * t + c[3]) * t + c[2]) * t + c[1]) / a + w
end


function alnrel(a::dp)::dp
    #-----------------------------------------------------------------------
    #            EVALUATION OF THE FUNCTION LN(1 + A)
    #-----------------------------------------------------------------------
    
    p = [-1.29418923021993, 0.405303492862024, -0.0178874546012214]
    q = [-16.2752256355323, 7.47811014037616, -0.0845104217945565]
    
    if abs(a) <= 0.375
        t = a / (a + 2.0)
        t2 = t * t
        w = (((p[3] * t2 + p[2]) * t2 + p[1]) * t2 + 1.0) / (((q[3] * t2 + q[2]) * t2 + q[1]) * t2 + 1.0)
        return 2.0 * t * w
    end
    
    x = 1.0 + a
    if a < 0.0
        x = (a + 0.5) + 0.5
    end
    return log(x)
end

function gamln1(a::dp)::dp
    #-----------------------------------------------------------------------
    #     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
    #-----------------------------------------------------------------------
    
    p = [0.577215664901533, 0.844203922187225, -0.168860593646662, -0.780427615533591,
         -0.402055799310489, -0.0673562214325671, -0.00271935708322958]
    q = [28.8743195473681, 31.2755088914843, 15.6875193295039, 3.61951990101499,
         0.325038868253937, 0.00667465618796164]
    r = [0.422784335098467, 0.848044614534529, 0.565221050691933, 0.156513060486551,
         0.017050248402265, 0.000497958207639485]
    s = [12.4313399877507, 5.48042109832463, 1.01552187439830, 0.0713309612391000,
         0.00116165475989616]
    
    if a < 0.6
        w = ((((((p[7] * a + p[6]) * a + p[5]) * a + p[4]) * a + p[3]) * a + p[2]) * a + p[1]) * a + p[1]
        return -a * w
    end
    
    x = (a - 0.5) - 0.5
    w = (((((r[6] * x + r[5]) * x + r[4]) * x + r[3]) * x + r[2]) * x + r[1]) * x + r[1]
    return x * w
end



end # module LogBeta

using .LogBeta

dp = Float64

function betain(x::dp, p::dp, q::dp, beta::dp)::dp
    #   Algorithm AS 63  Appl. Statist. (1973), vol.22, no.3
    #
    #   Computes incomplete beta function ratio for arguments
    #   x between zero and one, p and q positive.
    #   Log of complete beta function, beta, is assumed to be known
    #
    # ELF90-compatible version by Alan Miller
    # N.B. Argument IFAULT has been removed
    #
    # Latest revision - 5 July 2003

    zero = 0.0
    one = 1.0
    acu = 1.0e-14

    fn_val = x

    #   Test for admissibility of arguments

    if p <= zero || q <= zero
        println("AS63: Either p or q <= 0")
        return fn_val
    end
    if x < zero || x > one
        println("AS63: Argument x outside range (0, 1)")
        return fn_val
    end
    if x == zero || x == one
        return fn_val
    end

    #   Change tail if necessary and determine s

    psq = p + q
    cx = one - x
    if p < psq * x
        xx = cx
        cx = x
        pp = q
        qq = p
        indx = true
    else
        xx = x
        pp = p
        qq = q
        indx = false
    end
    term = one
    ai = one
    fn_val = one
    ns = qq + cx * psq

    #     Use Soper's reduction formulae.

    rx = xx / cx
    temp = qq - ai
    if ns == 0
        rx = xx
    end
    while true
        term = term * temp * rx / (pp + ai)
        fn_val += term
        temp = abs(term)
        if temp <= acu && temp <= acu * fn_val
            break
        end
        ai += one
        ns -= 1
        if ns >= 0
            continue
        end
        temp = psq
        psq += one
    end

    #     Calculate result

    fn_val = fn_val * exp(pp * log(xx) + (qq - one) * log(cx) - beta) / pp
    if indx
        fn_val = one - fn_val
    end

    return fn_val
end

using .LogBeta

a = 5.0
b = 3.0
log_beta_value = LogBeta.betaln(a, b)
println("The logarithm of the beta function for a = $a and b = $b is: $log_beta_value")

x = 0.5
p = 2.0
q = 3.0
beta = LogBeta.betaln(p, q)
incomplete_beta_value = betain(x, p, q, beta)
println("The incomplete beta function ratio for x = $x, p = $p, and q = $q is: $incomplete_beta_value")
