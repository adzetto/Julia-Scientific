module DipTest

function diptst(x::Vector{Float64}, n::Int, dip::Ref{Float64}, xl::Ref{Float64}, xu::Ref{Float64}, 
                ifault::Ref{Int}, gcm::Vector{Int}, lcm::Vector{Int}, mn::Vector{Int}, mj::Vector{Int})
    
    zero = 0.0
    half = 0.5
    one = 1.0

    ifault[] = 1
    if n <= 0
        return
    end
    ifault[] = 0

    # Check if N = 1
    if n != 1
        
        # Check that X is sorted
        ifault[] = 2
        for k in 2:n
            if x[k] < x[k-1]
                return
            end
        end
        ifault[] = 0
        
        # Check for all values of X identical,
        # and for 1 < N < 4.
        if x[n] > x[1] && n >= 4
            goto jump20
        end
    end
    xl[] = x[1]
    xu[] = x[n]
    dip[] = zero
    return

    # LOW contains the index of the current estimate of the lower end
    # of the modal interval, HIGH contains the index for the upper end.
    @label jump20
    fn = n
    low = 1
    high = n
    dip[] = one / fn
    xl[] = x[low]
    xu[] = x[high]

    # Establish the indices over which combination is necessary for the
    # convex minorant fit.
    mn[1] = 1
    for j in 2:n
        mn[j] = j - 1
        while true
            mnj = mn[j]
            mnmnj = mn[mnj]
            a = mnj - mnmnj
            b = j - mnj
            if mnj != 1 && (x[j] - x[mnj]) * a >= (x[mnj] - x[mnmnj]) * b
                mn[j] = mnmnj
            else
                break
            end
        end
    end

    # Establish the indices over which combination is necessary for the
    # concave majorant fit.
    mj[n] = n
    na = n - 1
    for jk in 1:na
        k = n - jk
        mj[k] = k + 1
        while true
            mjk = mj[k]
            mjmjk = mj[mjk]
            a = mjk - mjmjk
            b = k - mjk
            if mjk != n && (x[k] - x[mjk]) * a >= (x[mjk] - x[mjmjk]) * b
                mj[k] = mjmjk
            else
                break
            end
        end
    end

    # Start the cycling.
    # Collect the change points for the GCM from HIGH to LOW.
    @label jump70
    ic = 1
    gcm[1] = high

    while true
        igcm1 = gcm[ic]
        ic += 1
        gcm[ic] = mn[igcm1]
        if gcm[ic] <= low
            break
        end
    end
    icx = ic

    # Collect the change points for the LCM from LOW to HIGH.
    ic = 1
    lcm[1] = low
    while true
        lcm1 = lcm[ic]
        ic += 1
        lcm[ic] = mj[lcm1]
        if lcm[ic] >= high
            break
        end
    end
    icv = ic

    # ICX, IX, IG are counters for the convex minorant,
    # ICV, IV, IH are counters for the concave majorant.
    ig = icx
    ih = icv

    # Find the largest distance greater than 'DIP' between the GCM and
    # the LCM from LOW to HIGH.
    ix = icx - 1
    iv = 2
    d = zero
    if icx == 2 && icv == 2
        d = one / fn
        goto jump120
    end

    while true
        igcmx = gcm[ix]
        lcmiv = lcm[iv]
        if igcmx <= lcmiv
            
            # If the next point of either the GCM or LCM is from the LCM,
            # calculate the distance here.
            lcmiv1 = lcm[iv - 1]
            a = lcmiv - lcmiv1
            b = igcmx - lcmiv1 - 1
            dx = (x[igcmx] - x[lcmiv1] * a) / (fn * (x[lcmiv] - x[lcmiv1])) - b / fn
            ix -= 1
            if dx >= d
                d = dx
                ig = ix + 1
                ih = iv
            end
        else
            
            # If the next point of either the GCM or LCM is from the GCM,
            # calculate the distance here.
            lcmiv = lcm[iv]
            igcm = gcm[ix]
            igcm1 = gcm[ix + 1]
            a = lcmiv - igcm1 + 1
            b = igcm - igcm1
            dx = a / fn - ((x[lcmiv] - x[igcm1]) * b) / (fn * (x[igcm] - x[igcm1]))
            iv += 1
            if dx >= d
                d = dx
                ig = ix + 1
                ih = iv - 1
            end
        end

        ix = max(ix, 1)
        iv = min(iv, icv)
        if gcm[ix] == lcm[iv]
            break
        end
    end

    @label jump120
    if d >= dip[]
        
        # Calculate the DIPs for the current LOW and HIGH.
        # The DIP for the convex minorant.
        dl = zero
        if ig != icx
            icxa = icx - 1
            for j in ig:icxa
                temp = one / fn
                jb = gcm[j + 1]
                je = gcm[j]
                if je - jb > 1 && x[je] != x[jb]
                    a = je - jb
                    const = a / (fn * (x[je] - x[jb]))
                    for jr in jb:je
                        b = jr - jb + 1
                        t = b / fn - (x[jr] - x[jb]) * const
                        temp = max(temp, t)
                    end
                end
                dl = max(dl, temp)
            end
        end

        # The DIP for the concave majorant.
        du = zero
        if ih != icv
            icva = icv - 1
            for k in ih:icva
                temp = one / fn
                kb = lcm[k]
                ke = lcm[k + 1]
                if ke - kb > 1 && x[ke] != x[kb]
                    a = ke - kb
                    const = a / (fn * (x[ke] - x[kb]))
                    for kr in kb:ke
                        b = kr - kb - 1
                        t = (x[kr] - x[kb]) * const - b / fn
                        temp = max(temp, t)
                    end
                end
                du = max(du, temp)
            end
        end

        # Determine the current maximum.
        dipnew = max(dl, du)
        if dip[] < dipnew
            dip[] = dipnew
            low = gcm[ig]
            high = lcm[ih]
        end

        # Recycle
        goto jump70
    end

    dip[] = half * dip[]
    xl[] = x[low]
    xu[] = x[high]
end

end  # module DipTest
