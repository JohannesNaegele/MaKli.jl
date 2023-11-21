using MaKli
using Gadfly
using Distributions
import Plots

#--- Notiz ---#
# Bei der Aufgabe 1 kommt Quatsch raus.
# Zum Beispiel bei S(ω) für Ornstein-Uhlenbeck:
# julia> sol
# 100-element Vector{Float64}:
#  8.828580082817338
#  0.0245899118880504
#  0.000762129884829719
#  8.026327781434698e-5
#  1.082768773207192e-5
#  1.3696665677368428e-6
#  1.3840635450045305e-7
#  ⋮
#  0.0
#  0.0
# Wir haben die richtigen Parameter bei der Aufgabe 2 nicht gefunden.

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

T = 10
n = 100
N = 100
indices = 1:N
ω = [(i - 1) / N * π for i in indices]
sol = zeros(length(ω))
trajectory = zeros(n + 1)
function monte_carlo!(sol, trajectory, ω, f, g)
    sol .= 0.0
    temp_sol = zeros(ComplexF64, length(sol))
    for _ in 1:n^2
        temp_sol .= 0.0
        sde_euler!(trajectory, (u, t) -> f(u), (u, t) -> g(u), 1.0, (0.0, T), T / n)
        for i in eachindex(temp_sol)
            for k in eachindex(trajectory)
                temp_sol[i] += trajectory[k] * exp(-1im * π * (k-1) * ω[i])
                # temp_sol[i] += trajectory[k] * exp(-(i - 1) * π * (k-1) * ω[i])
            end
        end
        sol .+= abs.(temp_sol).^2 ./ n^2
        # sol .+= temp_sol.^2 ./ n^2
    end
    sol .*= T/(2π * n^2)
end

monte_carlo!(sol, trajectory, ω, ou_f, ou_g)

cind = [i for i in 0:(length(ω)-1)]
Plots.plot(cind[2:end], sol[2:end], xaxis=:log, yaxis=:log);
Plots.plot!(cind[2:end], (1 ./ (ω.^2))[2:end], xaxis=:log, yaxis=:log)

monte_carlo!(sol, trajectory,  ω, gb_f, gb_g)
Plots.plot(cind, sol, xaxis=:log, yaxis=:log);
Plots.plot!(cind, 1 ./ (1 .+ ω.^2), xaxis=:log, yaxis=:log)

# Aufgabe 2
const celsius_in_kelvin = 273.15
const S₀ = 1361 # Solar constant
const σ = 5.670374419e-8 # Stefan-Boltzmann constant
const H = 8007 # Skalenhöhe
const ρ = 1.2 # Dichte
const C = 1000
const ω_S = 2π / (60 * 60 * 24 * 365)

a_upper = 0.25
a_lower = 0.6
T_upper = 290
T_lower = 250

Q(t) = S₀ * (1 - 0.2 * sin(ω_S * t)) / (4H * ρ * C)

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

ε = 0.62 # Emissionsgrad
b = ε * σ / (H * ρ * C)
time = 3600 * 24 * 365 * 5 # 5 Jahr
n = 1e6
dt = time/1e6
trange = 0:dt:time
sol = zeros(length(trange))
sde_euler!(sol, u, (T, t) -> 0.0001, celsius_in_kelvin, (0.0, time), dt)
Plots.plot(sol)
sde_euler!(sol, u, (T, t) -> 0.001, celsius_in_kelvin, (0.0, time), dt)
Plots.plot(sol)
sde_euler!(sol, u, (T, t) -> 0.008, celsius_in_kelvin, (0.0, time), dt)
Plots.plot(sol)