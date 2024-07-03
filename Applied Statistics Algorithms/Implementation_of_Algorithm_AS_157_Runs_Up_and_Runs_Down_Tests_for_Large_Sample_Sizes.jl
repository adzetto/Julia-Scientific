module RunsTest

function udruns(x::Vector{Float64}, n::Int)
    ifault = 0
    uv = 0.0
    dv = 0.0

    # Define matrices
    a = [
        4529.4  9044.9 13568.0 18091.0 22615.0 27892.0;
        9044.9 18097.0 27139.0 36187.0 45234.0 55789.0;
        13568.0 27139.0 40721.0 54281.0 67852.0 83685.0;
        18091.0 36187.0 54281.0 72414.0 90470.0 111580.0;
        22615.0 45234.0 67852.0 90470.0 113262.0 139476.0;
        27892.0 55789.0 83685.0 111580.0 139476.0 172860.0
    ]
    b = [1.0 / 6.0, 5.0 / 24.0, 11.0 / 120.0, 19.0 / 720.0, 29.0 / 5040.0, 1.0 / 840.0]
    ucount = fill(0, 6)
    dcount = fill(0, 6)

    # Check the sample size
    if n < 4000
        ifault = n
        return uv, dv, ifault
    end

    ru = 1
    rd = 1
    for j in 2:n
        if x[j] < x[j-1]
            ucount[ru] += 1
            ru = 1
            if rd < 6 rd += 1 end
        elseif x[j] == x[j-1]
            ifault = 1
            return uv, dv, ifault
        else
            dcount[rd] += 1
            rd = 1
            if ru < 6 ru += 1 end
        end
    end
    ucount[ru] += 1
    dcount[rd] += 1

    # Calculate UV and DV
    for i in 1:6
        for j in 1:6
            uv += (ucount[i] - n * b[i]) * (ucount[j] - n * b[j]) * a[i, j]
            dv += (dcount[i] - n * b[i]) * (dcount[j] - n * b[j]) * a[i, j]
        end
    end
    uv /= n
    dv /= n

    return uv, dv, ifault
end

end # module
