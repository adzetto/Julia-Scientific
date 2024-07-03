module KMeans

export kmns

# Define constants
const zero = 0.0
const one = 1.0
const big = typemax(Float64)

function kmns(a::Matrix{Float64}, m::Int, n::Int, c::Matrix{Float64}, k::Int,
              ic1::Vector{Int}, ic2::Vector{Int}, nc::Vector{Int}, an1::Vector{Float64},
              an2::Vector{Float64}, ncp::Vector{Int}, d::Vector{Float64},
              itran::Vector{Int}, live::Vector{Int}, iter::Int, wss::Vector{Float64})

    ifault = 3
    if k <= 1 || k >= m
        return ifault
    end

    # For each point i, find its two closest centers, ic1(i) and ic2(i). Assign it to ic1(i).
    for i in 1:m
        ic1[i] = 1
        ic2[i] = 2
        dt = [zero, zero]
        for il in 1:2
            for j in 1:n
                da = a[i, j] - c[il, j]
                dt[il] += da * da
            end
        end
        if dt[1] > dt[2]
            ic1[i] = 2
            ic2[i] = 1
            temp = dt[1]
            dt[1] = dt[2]
            dt[2] = temp
        end

        for l in 3:k
            db = zero
            for j in 1:n
                dc = a[i, j] - c[l, j]
                db += dc * dc
                if db >= dt[2]
                    break
                end
            end
            if db >= dt[1]
                dt[2] = db
                ic2[i] = l
            else
                dt[2] = dt[1]
                ic2[i] = ic1[i]
                dt[1] = db
                ic1[i] = l
            end
        end
    end

    # Update cluster centers to be the average of points contained within them.
    for l in 1:k
        nc[l] = 0
        c[l, :] .= zero
    end

    for i in 1:m
        l = ic1[i]
        nc[l] += 1
        c[l, :] .+= a[i, :]
    end

    # Check if there is any empty cluster at this stage
    for l in 1:k
        if nc[l] == 0
            ifault = 1
            return ifault
        end
        aa = nc[l]
        c[l, :] ./= aa

        # Initialize AN1, AN2, ITRAN & NCP
        an2[l] = aa / (aa + one)
        an1[l] = big
        if aa > one
            an1[l] = aa / (aa - one)
        end
        itran[l] = 1
        ncp[l] = -1
    end

    indx = 0
    for ij in 1:iter
        # In this stage, there is only one pass through the data.
        # Each point is re-allocated, if necessary, to the cluster that will
        # induce the maximum reduction in within-cluster sum of squares.
        indx = optra!(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, indx)

        # Stop if no transfer took place in the last M optimal transfer steps.
        if indx == m
            break
        end

        # Each point is tested in turn to see if it should be re-allocated
        # to the cluster to which it is most likely to be transferred, IC2(I), 
        # from its present cluster, IC1(I). Loop through the data until no further change is to take place.
        qtran!(a, m, n, c, ic1, ic2, nc, an1, an2, ncp, d, itran, indx)

        # If there are only two clusters, there is no need to re-enter the optimal transfer stage.
        if k == 2
            break
        end

        # NCP has to be set to 0 before entering OPTRA.
        fill!(ncp, 0)
    end

    # Compute within-cluster sum of squares for each cluster.
    fill!(wss, zero)
    for l in 1:k
        c[l, :] .= zero
    end
    for i in 1:m
        ii = ic1[i]
        c[ii, :] .+= a[i, :]
    end
    for j in 1:n
        for l in 1:k
            c[l, j] /= nc[l]
        end
        for i in 1:m
            ii = ic1[i]
            da = a[i, j] - c[ii, j]
            wss[ii] += da * da
        end
    end

    ifault = 0
    return ifault
end

function optra!(a::Matrix{Float64}, m::Int, n::Int, c::Matrix{Float64}, k::Int,
                ic1::Vector{Int}, ic2::Vector{Int}, nc::Vector{Int}, an1::Vector{Float64},
                an2::Vector{Float64}, ncp::Vector{Int}, d::Vector{Float64},
                itran::Vector{Int}, live::Vector{Int}, indx::Int)
    
    for l in 1:k
        if itran[l] == 1
            live[l] = m + 1
        end
    end
    for i in 1:m
        indx += 1
        l1 = ic1[i]
        l2 = ic2[i]
        ll = l2

        # If point I is the only member of cluster L1, no transfer.
        if nc[l1] == 1
            if indx == m
                return indx
            end
            continue
        end

        # If L1 has not yet been updated in this stage, no need to re-compute D(I).
        if ncp[l1] == 0
            de = zero
            for j in 1:n
                df = a[i, j] - c[l1, j]
                de += df * df
            end
            d[i] = de * an1[l1]
        end

        # Find the cluster with minimum R2.
        da = zero
        for j in 1:n
            db = a[i, j] - c[l2, j]
            da += db * db
        end
        r2 = da * an2[l2]
        for l in 1:k
            # If I >= LIVE(L1), then L1 is not in the live set. 
            # If this is true, we only need to consider clusters that are in the live set for possible transfer of point I. 
            # Otherwise, we need to consider all possible clusters.
            if i >= live[l1] && i >= live[l] || l == l1 || l == ll
                continue
            end
            rr = r2 / an2[l]
            dc = zero
            for j in 1:n
                dd = a[i, j] - c[l, j]
                dc += dd * dd
                if dc >= rr
                    break
                end
            end
            if dc * an2[l] < r2
                r2 = dc * an2[l]
                l2 = l
            end
        end

        if r2 >= d[i]
            # If no transfer is necessary, L2 is the new IC2(I).
            ic2[i] = l2
            if indx == m
                return indx
            end
        else
            indx = 0
            live[l1] = m + i
            live[l2] = m + i
            ncp[l1] = i
            ncp[l2] = i
            al1 = nc[l1]
            alw = al1 - one
            al2 = nc[l2]
            alt = al2 + one
            for j in 1:n
                c[l1, j] = (c[l1, j] * al1 - a[i, j]) / alw
                c[l2, j] = (c[l2, j] * al2 + a[i, j]) / alt
            end
            nc[l1] -= 1
            nc[l2] += 1
            an2[l1] = alw / al1
            an1[l1] = big
            if alw > one
                an1[l1] = alw / (alw - one)
            end
            an1[l2] = alt / al2
            an2[l2] = alt / (alt + one)
            ic1[i] = l2
            ic2[i] = l1
        end
    end

    for l in 1:k
        itran[l] = 0
        live[l] -= m
    end

    return indx
end

function qtran!(a::Matrix{Float64}, m::Int, n::Int, c::Matrix{Float64}, ic1::Vector{Int},
                ic2::Vector{Int}, nc::Vector{Int}, an1::Vector{Float64}, an2::Vector{Float64},
                ncp::Vector{Int}, d::Vector{Float64}, itran::Vector{Int}, indx::Int)
    
    icoun = 0
    istep = 0
    while true
        for i in 1:m
            icoun += 1
            istep += 1
            l1 = ic1[i]
            l2 = ic2[i]

            # If point I is the only member of cluster L1, no transfer.
            if nc[l1] == 1
                continue
            end

            # If ISTEP > NCP(L1), no need to re-compute distance from point I to cluster L1.
            if istep > ncp[l1]
                da = zero
                for j in 1:n
                    db = a[i, j] - c[l1, j]
                    da += db * db
                end
                d[i] = da * an1[l1]
            end

            # If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of point I at this step.
            if istep >= ncp[l1] && istep >= ncp[l2]
                continue
            end

            r2 = d[i] / an2[l2]
            dd = zero
            for j in 1:n
                de = a[i, j] - c[l2, j]
                dd += de * de
                if dd >= r2
                    break
                end
            end

            if dd * an2[l2] < r2
                # Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters L1 & L2.
                # Also update IC1(I) & IC2(I). Note that if any updating occurs in this stage, INDX is set back to 0.
                icoun = 0
                indx = 0
                itran[l1] = 1
                itran[l2] = 1
                ncp[l1] = istep + m
                ncp[l2] = istep + m
                al1 = nc[l1]
                alw = al1 - one
                al2 = nc[l2]
                alt = al2 + one
                for j in 1:n
                    c[l1, j] = (c[l1, j] * al1 - a[i, j]) / alw
                    c[l2, j] = (c[l2, j] * al2 + a[i, j]) / alt
                end
                nc[l1] -= 1
                nc[l2] += 1
                an2[l1] = alw / al1
                an1[l1] = big
                if alw > one
                    an1[l1] = alw / (alw - one)
                end
                an1[l2] = alt / al2
                an2[l2] = alt / (alt + one)
                ic1[i] = l2
                ic2[i] = l1

                # If no re-allocation took place in the last M steps, return.
                if icoun == m
                    return
                end
            end
        end
    end
end

end # module KMeans

using .KMeans

# Example usage
m = 10
n = 2
k = 3
a = rand(m, n)
c = rand(k, n)
ic1 = zeros(Int, m)
ic2 = zeros(Int, m)
nc = zeros(Int, k)
an1 = zeros(Float64, k)
an2 = zeros(Float64, k)
ncp = zeros(Int, k)
d = zeros(Float64, m)
itran = zeros(Int, k)
live = zeros(Int, k)
iter = 100
wss = zeros(Float64, k)

ifault = KMeans.kmns(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, iter, wss)

println("Clusters: ", c)
println("Cluster assignments: ", ic1)
println("Within-cluster sum of squares: ", wss)
println("Ifault: ", ifault)
