using MaKli
using DifferentialEquations
using Plots
using Random

function sde_euler(f, g, u₀, tspan, dt)
    sol = Vector(undef, 1:dt:tspan[end])
    for i in eachindex(sol)[2:end]
        t = ...
        sol[i] = sol[i - 1] + f(sol[i - 1], t) * dt + g(sol[i - 1], t) * randn()
    end
    return sol
end

α = 2.0
β = 2.0
u₀ = 1.0
T = 200.0
f(u, p, t) = α * u
g(u, p, t) = β * u
dt = T / 10^3
tspan = (0.0, T)
prob = SDEProblem(f, g, u₀, tspan)
sol = solve(prob, dt = dt)
plot(sol)

function comparison(T)
    p = plot()
    for i in eachindex(α)
        Random.seed!(1)
    end

end