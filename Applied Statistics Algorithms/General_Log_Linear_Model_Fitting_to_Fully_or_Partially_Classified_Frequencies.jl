module GLLM

# Function to fit a generalized log-linear model
# to fully or partially classified frequencies.

function gllm(ni::Int, nid::Int, nj::Int, nk::Int, nkp::Int, ji::Vector{Int}, y::Vector{Float64}, 
              c::Array{Float64, 2}, conv::Float64, e::Vector{Float64}, f::Vector{Float64}, 
              cspr::Ref{Float64}, cslr::Ref{Float64}, ifault::Ref{Int})

    # Constants
    eps = 1e-5
    zero = 0.0
    one = 1.0

    # Initialize variables
    w = zeros(Float64, 4, ni)
    v = zeros(Float64, 2, nkp)
    cmin = zero
    ifault[] = 1

    # Check the ji array
    for i in 1:ni
        if ji[i] < 1 || ji[i] > nj
            return
        end
    end

    ifault[] = 0

    # Initialize
    e .= one
    w[3, 1:nj] .= zero

    # Standardize the C matrix
    for i in 1:ni
        for k in 1:nk
            cmin = min(cmin, c[i, k])
        end
    end

    if cmin != zero
        c .= c .- cmin
    end

    while true
        ctmax = zero
        for i in 1:ni
            csum = zero
            for k in 1:nk
                csum += c[i, k]
            end
            ctmax = max(ctmax, csum)
            w[4, i] = csum
        end

        if ctmax <= eps
            ifault[] = 2
            return
        end

        if abs(ctmax - one) > eps
            for i in 1:ni
                for k in 1:nk
                    c[i, k] /= ctmax
                end
            end
        else
            break
        end
    end

    for i in 1:ni
        if abs(w[4, i] - one) > eps
            goto nkk_not_equal_nk
        end
    end

    nkk = nk
    goto enter_em_algorithm

    nkk_not_equal_nk:
    nkk = nkp
    for i in 1:ni
        c[i, nkk] = one - w[4, i]
    end

    enter_em_algorithm:
    f .= zero
    for i in 1:ni
        j = ji[i]
        f[j] += e[i]
    end

    # Check for convergence
    lconv = true
    for j in 1:nj
        if abs(f[j] - w[3, j]) > conv
            lconv = false
        end
        w[3, j] = f[j]
    end

    if !lconv
        for i in 1:ni
            j = ji[i]
            w[1, i] = y[j]
            if f[j] > eps
                w[1, i] = e[i] * y[j] / f[j]
            end
        end

        for k in 1:nkk
            v[1, k] = sum(c[:, k] .* w[1, :])
        end

        while true
            w[2, :] .= e
            for k in 1:nkk
                v[2, k] = sum(c[:, k] .* e)
                for i in 1:ni
                    if c[i, k] > eps && v[2, k] > eps
                        e[i] *= (v[1, k] / v[2, k])^c[i, k]
                    end
                end
            end

            converged = true
            for i in 1:ni
                if abs(e[i] - w[2, i]) > conv
                    converged = false
                    break
                end
            end

            if converged
                break
            end
        end

        goto enter_em_algorithm
    end

    # Calculate the goodness-of-fit statistics
    cspr[] = zero
    cslr[] = zero
    for j in 1:nj
        yy = y[j]
        ff = f[j]
        if ff > eps
            cspr[] += (yy - ff)^2 / ff
            if yy > eps
                cslr[] += yy * log(yy / ff)
            end
        end
    end

    cslr[] *= 2.0
end

end  # module GLLM
