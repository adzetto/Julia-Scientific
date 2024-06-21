# Filename: random_force_particle.jl
# Description: Perform a 1000-second simulation of a particle's motion in the \(xy\)-plane while it is subjected to a random force.
# Author: adzetto
# Date: 2024-06-21

#=
16. Write a program to simulate the following phenomenon: A particle is moving in the $x y$-plane under the effect of a random force. It starts at $(0,0)$. At the end of each second, it moves 1 unit in a random direction. We want to record in a table its position at the end of each second, taking altogether 1000 seconds.
=#

using Random
using Printf

function random_step()
    angle = 2 * π * rand()
    return cos(angle), sin(angle)
end

function simulate_particle_motion(steps::Int)
    x, y = 0.0, 0.0
    positions = Vector{Tuple{Float64, Float64}}()

    for _ in 1:steps
        dx, dy = random_step()
        x += dx
        y += dy
        push!(positions, (x, y))
    end

    return positions
end

num_steps = 1000

positions = simulate_particle_motion(num_steps)

println("Time\tX\t\tY")
for (i, (x, y)) in enumerate(positions)
    @printf("%d\t%.6f\t%.6f\n", i, x, y)
end
