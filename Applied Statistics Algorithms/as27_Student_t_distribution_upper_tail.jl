# N.B. Argument IFAULT has been removed.
# Date: 2024-07-02  Time: 10:12:00
# ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
# Calculate the upper tail area under Student's t-distribution
# Translated from FORTRAN by adzetto

function studnt(t, doff)
    # Local variables
    v, x, tt = 0.0, 0.0, 0.0
    pos = false
    four, one, zero, half = 4.0, 1.0, 0.0, 0.5

    a1, a2, a3, a4, a5 = 0.09979441, -0.581821, 1.390993, -1.222452, 2.151185
    b1, b2 = 5.537409, 11.42343
    c1, c2, c3, c4, c5 = 0.04431742, -0.2206018, -0.03317253, 5.679969, -12.96519
    d1, d2 = 5.166733, 13.49862
    e1, e2, e3, e4, e5 = 0.009694901, -0.1408854, 1.88993, -12.75532, 25.77532
    f1, f2 = 4.233736, 14.3963
    g1, g2, g3, g4, g5 = -9.187228e-5, 0.03789901, -1.280346, 9.249528, -19.08115
    h1, h2 = 2.777816, 16.46132
    i1, i2, i3, i4, i5 = 5.79602e-4, -0.02763334, 0.4517029, -2.657697, 5.127212
    j1, j2 = 0.5657187, 21.83269

    # Check that number of degrees of freedom > 4.
    if doff <= four
        println("** Error in AS27 - degrees of freedom <= 4  **")
        return
    end

    # Evaluate series.
    v = one / doff
    pos = (t >= zero)
    tt = abs(t)
    x = half * (one + 
        tt * ((a1 + v * (a2 + v * (a3 + v * (a4 + v * a5)))) / (one - v * (b1 - v * b2)) + 
        tt * ((c1 + v * (c2 + v * (c3 + v * (c4 + v * c5)))) / (one - v * (d1 - v * d2)) + 
        tt * ((e1 + v * (e2 + v * (e3 + v * (e4 + v * e5)))) / (one - v * (f1 - v * f2)) + 
        tt * ((g1 + v * (g2 + v * (g3 + v * (g4 + v * g5)))) / (one - v * (h1 - v * h2)) + 
        tt * ((i1 + v * (i2 + v * (i3 + v * (i4 + v * i5)))) / (one - v * (j1 - v * j2))))))))^(-8)

    if pos
        return x
    else
        return one - x
    end
end


t = 2.0
doff = 5.0

# Calculate the upper tail area under Student's t-distribution
result = studnt(t, doff)

println("The upper tail area under Student's t-distribution for t = $t and degrees of freedom = $doff is $result")
