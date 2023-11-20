using MaKli
using Gadfly

function sde_euler!(sol, f, g, u₀, tspan, dt)
    time = 0:dt:tspan[end]
    sol[1] = u₀
    for t in eachindex(time)[2:end]
        t_value = time[t]
        increment = rand(Normal(0, sqrt(dt)))
        sol[t] = sol[t - 1] + f(sol[t - 1], t_value) * dt + g(sol[t - 1], t_value) * increment
    end
    return time, sol
end

# Aufgabe 1

ou_f(u) = 0.1
ou_g(u) = 1.0
gb_f(u) = 0.1
gb_g(u) = 0.2

sde_euler!(sol, (u, t) -> ou_f(u), (u, t) -> ou_g(u), 1.0, (0.0, T), T / n)

# Aufgabe 2
const celsius_in_kelvin = 273.15
const S₀ = 1361 # Solar constant
const ε = 0.62 # Emissionsgrad
const σ = 5.670374419e-8 # Stefan-Boltzmann constant
const H = 8007 # Skalenhöhe
const ρ = 1.2 # Dichte
const C = 1000
const ω_S = 2π / (3600 * 24 * 365)
const b = ε * σ / (H * ρ * C)

a_upper = 0.25
a_lower = 0.6
T_upper = 290
T_lower = 260

Q(t) = S₀ * (1 - 0.06 * sin(ω_S * t)) / (4H * ρ * C)

function α(T)
    if T > T_upper
        return a_upper
    elseif T < T_lower
        return a_lower
    else
        return (T - T_lower) / (T_upper - T_lower) * (a_upper - a_lower) + a_lower
    end
end

R_in(T, Q) = Q * (1 - α(T))
R_out(T) = b * T^4
u(T, t) = R_in(T, Q(t)) - R_out(T)

# anim = Plots.@animate for t in [(1 / ω_S) * (i / 100) * 2π for i in 1:100]
#     plot(temperatures, map(T -> R_in(T, Q(t)), temperatures), dpi=200, label="R_in", xlabel="Temperature in Kelvin", ylabel="Radiation in W/m²")
#     plot!(temperatures, map(T -> R_out(T), temperatures), label="R_out")
# end

sde_euler!(sol, (T, t) -> u(T, t), (T, t) -> ε, 1.0, (0.0, T), T / n)
