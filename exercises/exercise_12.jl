using MaKli
using Gadfly
using Distributions
import Plots

# Aufgabe 1

# (a)

# (b)

# (c)

# Aufgabe 2

function sde_euler!(f, g, u₀, tspan, dt)
    t0, t_end = tspan
    time = t0:dt:t_end
    sol = zeros(length(u₀), length(time))
    # increment = zeros(first(methods(g, typeof.(u₀))).nargs - 1)
    increment = similar(g(u₀, t0))
    sol[:, 1] = u₀
    for t in eachindex(time)[2:end]
        t_value = time[t]
        for i in eachindex(increment)
            if increment[i] isa Number
                increment[i] = rand(Normal(0, sqrt(dt)))
            else
                increment[i] = rand(Normal(0, sqrt(dt)), length(increment[i]))
            end
        end
        sol[:, t] = sol[:, t - 1] .+ f(sol[:, t - 1], t_value) .* dt .+ g(sol[:, t - 1], t_value) .* increment
    end
    return time, sol
end

const ε = 0.1
const T = 100
n = 1e7

# (a)
f(x, y; a, b) = [
    y - x,
    -a * x + b * y - x^2 * y
]

g(x, y) = [
    0.0,
    ε
]

res = sde_euler!((u, t) -> f(u[1], u[2], a=1.4, b=3.0), (u, t) -> g(u[1], u[2]), [1.0, 1.0], (0.0, T), T / n)
Plots.plot(res[2][1, :], res[2][2, :]);
Plots.savefig("./exercises/graphics/12_2a_11.png")
Plots.plot(res[1], res[2][1, :], res[2][2, :]);
Plots.savefig("./exercises/graphics/12_2a_12.png")

# (b)
res = sde_euler!((u, t) -> f(u[1], u[2], a=4.4, b=0.3), (u, t) -> g(u[1], u[2]), [1.0, 1.0], (0.0, T), T / n)
Plots.plot(res[2][1, :], res[2][2, :]);
Plots.savefig("./exercises/graphics/12_2a_21.png")
Plots.plot(res[1], res[2][1, :], res[2][2, :]);
Plots.savefig("./exercises/graphics/12_2a_22.png")