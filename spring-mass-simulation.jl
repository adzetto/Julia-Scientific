using Plots
using LaTeXStrings
using Printf

m = 30.0
k = 50.0
mu_s = 0.05
mu_d = 0.05
x0 = 6.0
v0 = 0.0
T = 25.0
dt = 1e-4

g = 9.81
n = ceil(Int, T / dt) + 1
t = range(0, stop = T, length = n)

x = zeros(n)
v = zeros(n)

x[1] = x0
v[1] = v0

function fr(x, v, m, g, k, mu_s, mu_d)
    if abs(v) > 1e-20
        return mu_d * m * g * sign(v)
    else
        return -min(mu_s * m * g, abs(k * x)) * sign(x)
    end
end

for i in 2:n
    x[i] = x[i-1] + dt * v[i-1]
    v[i] = v[i-1] + dt * (-1 / m * fr(x[i-1], v[i-1], m, g, k, mu_s, mu_d) - k / m * x[i-1])
end

p1 = plot(t, x, label=L"x(t)", xlabel=L"Time\, t\, [s]", ylabel=L"Displacement\, x\, [m]", title="Displacement vs. Time", legend=:bottomright)
p2 = plot(x, v, label=L"v(t)", xlabel=L"Displacement\, x\, [m]", ylabel=L"Velocity\, v\, [m/s]", title="Velocity vs. Displacement", legend=:bottomright)

plot(p1, p2, layout=(1, 2), size=(1000, 400), titlefontsize=10, guidefontsize=10, tickfontsize=8)
savefig("1.pdf")

m = 30.0
k = 50.0
g = 9.81
mu = 0.05
x0_initial = 6.0
N_cyc = 5

function solve_quadratic(A, B, C)
    discriminant = B^2 - 4 * A * C
    if discriminant < 0
        error("Discriminant is negative. No real roots.")
    end
    root1 = (-B + sqrt(discriminant)) / (2 * A)
    root2 = (-B - sqrt(discriminant)) / (2 * A)
    return root1, root2
end

results = []

x0 = x0_initial
for cycle in 1:N_cyc
    A, B = 0.5 * k, m * g * mu
    C0 = m * g * mu * x0 - 0.5 * k * x0^2
    x1_1, x1_2 = solve_quadratic(A, B, C0)
    x1 = 0 < x1_1 < x0 ? x1_1 : x1_2

    C1 = m * g * mu * x1 - 0.5 * k * x1^2
    x2_1, x2_2 = solve_quadratic(A, B, C1)
    x2 = 0 < x2_1 < x1 ? x2_1 : x2_2

    push!(results, (cycle - 0.5, abs(x0), abs(x1)))
    push!(results, (cycle, abs(x1), abs(x2)))

    x0 = x2
end

for i in eachindex(results)
    if i % 2 == 0
        results[i] = (results[i][1], results[i][2], -results[i][3])
    else
        results[i] = (results[i][1], -results[i][2], results[i][3])
    end
end

cycle_numbers = []
positions = []

for (cycle, start_pos, end_pos) in results
    push!(cycle_numbers, cycle)
    push!(positions, start_pos)
    push!(cycle_numbers, cycle)
    push!(positions, end_pos)
end

colors = cgrad(:viridis, length(cycle_numbers))

fig = plot(size=(800, 600), title="Position of Mass Over Cycles", xlabel="Cycle", ylabel="Position [m]", legend=false, grid=:both)
for i in 1:2:length(cycle_numbers)-1
    plot!(cycle_numbers[i:i+1], positions[i:i+1], marker=:circle, color=colors[i], linewidth=2)
end
for i in 1:2:length(positions)
    annotate!([cycle_numbers[i]], [positions[i]], text(@sprintf("%.4f", positions[i]), 10, :center))
end
cbar = Colorbar(fig, cgrad(:viridis, length(cycle_numbers)), ticks=1:N_cyc, label="Cycle Number")
display(fig)
savefig(fig, "2.pdf")

function compute_velocity_profile(x0, x, m, k, g, mu)
    return sqrt(max(0, (-k * x^2 + 2 * (k * x0 - m * g * mu) * x) / m))
end

velocities = [compute_velocity_profile(x0_initial, x, m, k, g, mu) for x in positions]

p3 = plot(cycle_numbers, velocities, marker=:x, title="Velocity Profile of Mass Block Over Cycles", xlabel="Cycle", ylabel="Velocity [m/s]", grid=:both)
savefig(p3, "3.pdf")
