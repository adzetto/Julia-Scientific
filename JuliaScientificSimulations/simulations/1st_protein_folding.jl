#=
Filename: protein_folding.jl
Description: Use gradient descent to minimise the Lennard-Jones potential and predict the protein's structure.
Author: adzetto
Date: 2024-06-21
=#


#=
10. (Molecular conformation: Protein folding project) Forces that govern folding of amino acids into proteins are due to bonds between individual atoms and to weaker interactions between unbound atoms such as electrostatic and Van der Waals forces. The Van der Waals forces are modeled by the Lennard-Jones potential$$U(r)=\frac{1}{r^{12}}-\frac{2}{r^6}$$where $r$ is the distance between atoms.In the figure, the energy minimum is -1 and it is achieved at $r=1$. Explore this subject and the numerical methods used. One| approach is to predict the conformation of the proteins in finding the minimum potential energy of the total configuration of amino acids. For a cluster of atoms with positions $\left(x_1, y_1, z_1\right)$ to $\left(x_n, y_n, z_n\right)$, the objective function to be minimized is$$U=\sum_{i<j} \frac{1}{r_{i j}^{12}}-\frac{2}{r_{i j}^6}$$over all pairs of atoms. Here, $r_{i j}=\left[\left(x_i-x_j\right)^2+\left(y_i-y_j\right)^2+\left(z_i-z_j\right)^2\right]^2$ is the distance between atoms $i$ and $j$. This optimization problem finds the rectangular coordinates of the atoms. See Sauer [2006] for additional details.
=#

using Random
using LinearAlgebra

# Lennard-Jones p.f.
function lennard_jones_potential(r)
    r12 = r^12
    r6 = r^6
    return 1 / r12 - 2 / r6
end

function distance(atom1, atom2)
    return norm(atom1 .- atom2)
end

function total_potential_energy(positions)
    n = size(positions, 1)
    U = 0.0
    for i in 1:n-1
        for j in i+1:n
            r_ij = distance(positions[i, :], positions[j, :])
            U += lennard_jones_potential(r_ij)
        end
    end
    return U
end

function gradient(positions)
    n = size(positions, 1)
    grad = zeros(size(positions))
    for i in 1:n-1
        for j in i+1:n
            r_ij = distance(positions[i, :], positions[j, :])
            r_vec = positions[i, :] .- positions[j, :]
            r2 = r_ij^2
            r8 = r2^4
            r14 = r2^7
            coeff = (12 / r14 - 12 / (2 * r8))
            grad[i, :] += coeff * r_vec
            grad[j, :] -= coeff * r_vec
        end
    end
    return grad
end

function gradient_descent(positions, learning_rate=0.001, max_iters=10000, tolerance=1e-6)
    for iter in 1:max_iters
        grad = gradient(positions)
        positions -= learning_rate * grad
        current_potential = total_potential_energy(positions)
        if norm(grad) < tolerance
            println("Converged after $iter iterations.")
            break
        end
        if iter % 100 == 0
            println("Iteration $iter: Potential Energy = $current_potential")
        end
    end
    return positions
end

n_atoms = 5
positions = rand(n_atoms, 3)

optimized_positions = gradient_descent(positions)

println("Optimized Positions:")
println(optimized_positions)

final_potential_energy = total_potential_energy(optimized_positions)
println("Final Potential Energy: $final_potential_energy")
