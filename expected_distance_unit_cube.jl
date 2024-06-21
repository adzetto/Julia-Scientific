# Filename: expected_distance_unit_cube.jl
# Description: Applying Monte Carlo simulation, determine the predicted distance between two randomly selected points in a unit cube.
# Author: adzetto
# Date: 2024-06-21

#=
23. In the unit cube $\{(x, y, z): 0 \leqq x \leqq 1,0 \leqq y \leqq 1,0 \leqq z \leqq 1\}$, if two points are randomly chosen, then what is the expected distance between them?
=#

using Random

function random_point()
    return rand(), rand(), rand()
end

function euclidean_distance(p1::Tuple{Float64, Float64, Float64}, p2::Tuple{Float64, Float64, Float64})
    return sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2 + (p1[3] - p2[3])^2)
end

function estimate_expected_distance(num_simulations::Int)
    total_distance = 0.0

    for _ in 1:num_simulations
        point1 = random_point()
        point2 = random_point()
        total_distance += euclidean_distance(point1, point2)
    end

    return total_distance / num_simulations
end

num_simulations = 100000000

expected_distance = estimate_expected_distance(num_simulations)

println("Estimated expected distance between two randomly chosen points in the unit cube: $expected_distance")
