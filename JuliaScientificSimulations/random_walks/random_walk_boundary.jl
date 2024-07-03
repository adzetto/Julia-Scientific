# Filename: random_walk_boundary.jl
# Description: Try to predict the likelihood that the particle will move randomly over a square and see whether it lands on the bottom first.
# Author: adzetto
# Date: 2024-06-21

#=
18. (Another random walk) Consider the lattice points (points with integer coordinates) in the square $0 \leqq x \leqq 6,0 \leqq y \leqq 6$. A particle starts at the point $(4,4)$ and moves in the following way: At each step, it moves with equal probability to one of the four adjacent lattice points. What is the probability that when the particle first crosses the boundary of the square, it crosses the bottom side? Use Monte Carlo simulation.
=#

using Random

function simulate_walk()
    x, y = 4, 4
    
    while 0 <= x <= 6 && 0 <= y <= 6
        r = rand()
        if r < 0.25
            x += 1
        elseif r < 0.5
            x -= 1
        elseif r < 0.75
            y += 1
        else
            y -= 1
        end

        if x < 0 || x > 6 || y < 0 || y > 6
            return y < 0
        end
    end
end

function estimate_probability(num_simulations::Int)
    count = 0
    
    for _ in 1:num_simulations
        if simulate_walk()
            count += 1
        end
    end
    
    return count / num_simulations
end

num_simulations = 1000000

probability = estimate_probability(num_simulations)

println("Estimated probability that the particle first crosses the bottom side: $probability")
