using LinearAlgebra
using Plots

function f(x)
    return [x[1]^2 + x[2]^2 - 1; x[1]^2 - x[2]]
end

function J(x)
    return [2*x[1] 2*x[2]; 2*x[1] -1]
end

function newtons_method(f, J, x0, tol=1e-8, max_iter=100)
    x = x0
    errors = []
    
    for k in 1:max_iter
        fx = f(x)
        Jx = J(x)
        
        dx = - Jx \ fx
        x = x + dx
        
        error = norm(f(x))
        push!(errors, error)
        
        if error < tol
            println("Converged in $k iterations.")
            return x, errors
        end
    end
    
    println("Did not converge within the maximum number of iterations.")
    return x, errors
end

x0 = [0.5; 0.5]
tol = 1e-8
max_iter = 100

root, errors = newtons_method(f, J, x0, tol, max_iter)

println("Found root: $root")

p = plot(errors, yscale=:log10, xlabel="Iteration", ylabel="Error", title="Newton's Method Convergence")
savefig(p, "/mnt/c/Users/lenovo/Documents/GitHub/Julia-Scientific/newton_convergence.pdf")
