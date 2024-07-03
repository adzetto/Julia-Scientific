module ContingencyTable

using Printf

function enum(r::Int, c::Int, n::Vector{Int}, m::Vector{Int})
    # Constants
    maxt = 10
    nmax = 201
    zero = 0.0
    one = 1

    # Initialize variables
    table = zeros(Int, maxt, maxt)
    bound = zeros(Int, maxt, maxt)
    z = zeros(Int, maxt)
    reps = zeros(Int, maxt)
    mult = ones(Int, maxt)
    rept = falses(Bool, maxt)
    reptc = falses(Bool, maxt)
    prob = zeros(Float64, maxt, maxt)
    flogm = zeros(Float64, nmax)
    factlm = zeros(Float64, nmax)

    # Check input values
    ifault = 1
    if r > maxt || c > maxt || r <= 0 || c <= 0
        return ifault
    end

    ifault = 3
    ntotal = sum(n)
    ntot = sum(m)
    
    ifault = 2
    if ntot != ntotal
        return ifault
    end

    ifault = 4
    if ntotal >= nmax
        return ifault
    end

    ifault = 0

    # Initialize flogm and factlm arrays
    flogm[1] = zero
    factlm[1] = zero
    for k in 1:ntotal
        flogm[k+1] = log(k)
        factlm[k+1] = factlm[k] + flogm[k+1]
    end

    rm = r - 1
    cm = c - 1

    # Sort rows and columns into ascending order
    sort!(n)
    sort!(m)

    # Calculate multiplicities of rows and columns
    multc = one
    repsc = one
    reptc[1] = false
    for j in 2:c
        reptc[j] = m[j] == m[j-1]
        if !reptc[j]
            repsc = one
        else
            repsc += one
            multc *= repsc
        end
    end

    multr = one
    repsr = one
    rept[1] = false
    for i in 2:r
        rept[i] = n[i] == n[i-1]
        if !rept[i]
            repsr = one
        else
            repsr += one
            multr *= repsr
        end
    end

    # If column multiplicity exceeds row multiplicity, transpose the table
    if multr < multc
        maxrc = max(r, c)
        for ij in 1:maxrc
            keep = n[ij]
            n[ij] = m[ij]
            m[ij] = keep
        end
        keep = r
        r = c
        c = keep
        rm = r - 1
        cm = c - 1
        for i in 1:r
            rept[i] = reptc[i]
        end
        multr = multc
    end

    # Set up initial table
    mult[1] = multr
    reps[1] = one
    prob0 = -factlm[ntotal+1]
    for i in 1:r
        prob0 += factlm[n[i]+1]
    end
    for j in 1:c
        prob0 += factlm[m[j]+1]
    end

    for j in 1:c
        bound[1, j] = m[j]
    end

    for i in 1:r
        if i != 1
            prob0 = prob[i-1, c]
        end
        left = n[i]
        ieqim = rept[i]
        for j in 1:cm
            ij = min(left, bound[i, j])
            table[i, j] = ij
            if j == 1
                prob[i, j] = prob0 - factlm[ij+1]
            else
                prob[i, j] = prob[i, j-1] - factlm[ij+1]
            end
            left -= table[i, j]
            if i < r
                bound[i+1, j] = bound[i, j] - table[i, j]
            end
            if left == 0
                goto 160
            end
            if ieqim
                ieqim = table[i, j] == table[i-1, j]
            end
        end
        table[i, c] = left
        prob[i, c] = prob[i, cm] - factlm[left+1]
        if i < r
            bound[i+1, c] = bound[i, c] - left
        end
        goto 180

        160:
        jp = j + 1
        for jj in jp:c
            table[i, jj] = 0
            prob[i, jj] = prob[i, jj-1]
            bound[i+1, jj] = bound[i, jj]
        end

        180:
        if i != 1
            mult[i] = mult[i-1]
            reps[i] = one
            if ieqim
                reps[i] = reps[i-1] + one
                mult[i] /= reps[i]
            end
        end
    end

    # Call eval for table 1
    eval(table, r, c, n, m, prob[r, c], mult[r])

    # Commence enumeration of remaining tables
    # Start of main loop
    while true
        i = r
        while i > 0
            i -= 1
            j = cm
            left = table[i, c]
            rowbnd = bound[i, c]
            while true
                if table[i, j] <= 0 || left >= rowbnd
                    if j == 1
                        break
                    end
                    left += table[i, j]
                    rowbnd += bound[i, j]
                    j -= 1
                else
                    ij = table[i, j]
                    prob[i, j] += flogm[ij+1]
                    table[i, j] -= 1
                    bound[i+1, j] += 1
                    if reps[i] != one
                        reps[i] = one
                        mult[i] = mult[i-1]
                    end
                    ii = i
                    iip = ii + 1
                    iim = ii - 1
                    jnext = j + 1
                    left += 1
                    goto 320
                end
            end
        end

        # If I = 0 no more tables are possible
        if i == 0
            return ifault
        end

        320:
        if jnext != c
            for j in jnext:cm
                table[ii, j] = min(left, bound[ii, j])
                left -= table[ii, j]
                bound[iip, j] = bound[ii, j] - table[ii, j]
                ij = table[ii, j]
                if j == 1
                    prob[ii, j] = prob[iim, c] - factlm[ij+1]
                else
                    prob[ii, j] = prob[ii, j-1] - factlm[ij+1]
                end
                if left == 0
                    goto 340
                end
            end
        end
        table[ii, c] = left
        prob[ii, c] = prob[ii, cm] - factlm[left+1]
        bound[iip, c] = bound[ii, c] - left
        goto 360

        340:
        jp = j + 1
        for jj in jp:c
            table[ii, jj] = 0
            prob[ii, jj] = prob[ii, jj-1]
            bound[iip, jj] = bound[ii, jj]
        end

        360:
        reps[ii] = one
        if ii > 1
            mult[ii] = mult[iim]
        end
        ii += 1
        iip = ii + 1
        iim = ii - 1
        if !rept[ii]
            left = n[ii]
            jnext = 1
        else
            for j in 1:c
                if bound[ii, j] > table[iim, j]
                    goto 250
                end
                ij = table[iim, j]
                table[ii, j] = ij
                bound[iip, j] = bound[ii, j] - table[ii, j]
                if j == 1
                    prob[ii, j] = prob[iim, c] - factlm[ij+1]
                else
                    prob[ii, j] = prob[ii, j-1] - factlm[ij+1]
                end
            end
            reps[ii] = reps[iim] + one
            mult[ii] = mult[iim] / reps[ii]
            goto 230
        end

        230:
        ii += 1
        if ii == r
            goto 370
        end
        iip = ii + 1
        iim = ii - 1
        if !rept[ii]
            left = n[ii]
            jnext = 1
        else
            for j in 1:c
                if bound[ii, j] > table[iim, j]
                    goto 250
                end
                ij = table[iim, j]
                table[ii, j] = ij
                bound[iip, j] = bound[ii, j] - table[ii, j]
                if j == 1
                    prob[ii, j] = prob[iim, c] - factlm[ij+1]
                else
                    prob[ii, j] = prob[ii, j-1] - factlm[ij+1]
                end
            end
            reps[ii] = reps[iim] + one
            mult[ii] = mult[iim] / reps[ii]
            goto 230
        end

        370:
        if !rept[r]
            ij = bound[r, 1]
            table[r, 1] = ij
            prob[r, 1] = prob[rm, c] - factlm[ij+1]
            for j in 2:c
                ij = bound[r, j]
                table[r, j] = ij
                prob[r, j] = prob[r, j-1] - factlm[ij+1]
            end
            mult[r] = mult[rm]
        else
            for j in 1:c
                if bound[r, j] > table[rm, j]
                    goto 400
                end
                ij = bound[r, j]
                table[r, j] = ij
                if j == 1
                    prob[r, j] = prob[rm, c] - factlm[ij+1]
                else
                    prob[r, j] = prob[r, j-1] - factlm[ij+1]
                end
                if table[r, j] != table[rm, j]
                    goto 410
                end
            end
            reps[r] = reps[rm] + one
            mult[r] = mult[rm] / reps[r]
            goto 430

            400:
            i = rm
            goto 210

            410:
            jp = j + 1
            for jj in jp:c
                ij = bound[r, jj]
                table[r, jj] = ij
                prob[r, jj] = prob[r, jj-1] - factlm[ij+1]
            end
            mult[r] = mult[rm]
        end

        430:
        eval(table, r, c, n, m, prob[r, c], mult[r])

        # End of main loop
    end
end

function eval(table::Array{Int, 2}, r::Int, c::Int, n::Vector{Int}, m::Vector{Int}, prob::Float64, mult::Int)
    # This is a dummy eval function. Modify it to perform the desired operations on each table.
    println("Table:")
    for i in 1:r
        for j in 1:c
            @printf("%4d", table[i, j])
        end
        println()
    end
    @printf("Probability: %.6f, Multiplicity: %d\n", prob, mult)
end

# Example usage
function main()
    r = 3
    c = 3
    n = [5, 3, 2]
    m = [4, 3, 3]
    ifault = enum(r, c, n, m)
    if ifault != 0
        @printf("IFAULT = %d\n", ifault)
    end
end

main()

end # module ContingencyTable
