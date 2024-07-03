module GCount

# Function to generate all possible N-bit binary codes and apply a user-defined procedure
function gcount(n::Int, apply::Function)
    if n < 1
        return 1  # ifault = 1
    end
    
    # Initialize variables
    status = falses(n)
    tpoint = collect(2:n+1)
    
    # Apply the user function to the initial state
    apply(n, n, status)

    while true
        # Generate new code
        change = 0
        if status[1]
            status[1] = false
            change = tpoint[2]
        else
            status[1] = true
            change = 2
        end
        apply(n, 1, status)

        # Check if count exhausted
        if change > n
            return 0  # ifault = 0
        end

        if status[change]
            status[change] = false
            tpoint[change] = tpoint[change + 1]
        else
            status[change] = true
            tpoint[change] = change + 1
        end
        apply(n, change, status)
    end
end

# Example user-defined procedure to print the binary codes
function print_comb(n::Int, change::Int, status::Vector{Bool})
    println(status)
end

end  # module GCount

using .GCount

# Example usage of the gcount function
function main()
    n = 4
    ifault = GCount.gcount(n, GCount.print_comb)
    if ifault != 0
        println("Error: ifault = ", ifault)
    end
end

main()
