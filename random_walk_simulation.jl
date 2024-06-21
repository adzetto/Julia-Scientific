# Filename: random_walk_simulation.jl
# Description: Simulates a random walk to estimate the probability of being more than 20 units away from the starting point after 50 steps.
# Author: adzetto
# Date: 2024-06-21

#=
17. (A random walk) On a windy night, a drunkard begins walking at the origin of a two-dimensional coordinate system. His steps are 1 unit in length and are random in the following way: With probability $\frac{1}{6}$, he takes a step east; with probability $\frac{1}{4}$, he takes a step north; with probability $\frac{1}{4}$, he takes a step south; and with probability $\frac{1}{3}$, he takes a step west. What is the probability that after 50 steps, he will be more than 20 units distant from the origin? Write a program to simulate this problem.
=#

using Random

function simulate_walk(steps::Int)
    x, y = 0, 0
    
    for _ in 1:steps
        r = rand()
        if r < 1/6
            x += 1
        elseif r < 1/6 + 1/4
            y += 1
        elseif r < 1/6 + 1/4 + 1/4
            y -= 1
        else
            x -= 1
        end
    end
    
    return hypot(x, y)
end

function estimate_probability(steps::Int, threshold::Float64, num_simulations::Int)
    count = 0
    
    for _ in 1:num_simulations
        distance = simulate_walk(steps)
        if distance > threshold
            count += 1
        end
    end
    
    return count / num_simulations
end

steps = 50
threshold = 20.0
num_simulations = 10e12

probability = estimate_probability(steps, threshold, num_simulations)

println("Estimated probability that the drunkard is more than $threshold units away from the origin after $steps steps: $probability")
