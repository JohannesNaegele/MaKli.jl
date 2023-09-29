using Pkg;
Pkg.activate("./MaKli");
using MaKli
using DifferentialEquations
import Gadfly
using DataFrames
import Cairo, Fontconfig
using Plots

const celsius_in_kelvin = 273.15
const S₀ = 1361 # Solar constant
const ε = 0.62 # Emissionsgrad
const σ = 5.670374419e-8 # Stefan-Boltzmann constant
const H = 8007 # Skalenhöhe
const ρ = 1.2 # Dichte
const C = 1000
const ω_S = 2π / (3600 * 24 * 365)
const b = ε * σ / (H * ρ * C)

const a_upper = 0.25
const a_lower = 0.6
const T_upper = 290
const T_lower = 260

# (a) Typo in t
function α(T)
    if T > T_upper
        return a_upper
    elseif T < T_lower
        return a_lower
    else
        return (T - T_lower) / (T_upper - T_lower) * (a_upper - a_lower) + a_lower
    end
end
Q(t) = S₀ * (1 - 0.06 * sin(ω_S * t)) / (4H * ρ * C)
R_in(T, Q) = Q * (1 - α(T))
R_out(T) = b * T^4
u(T, t) = R_in(T, Q(t)) - R_out(T)

temperatures = (-30:0.01:30) .+ celsius_in_kelvin

gr(windowsize=(480, 480))
anim = Plots.@animate for t in [(1 / ω_S) * (i / 100) * 2π for i in 1:100]
    plot(temperatures, map(T -> R_in(T, Q(t)), temperatures), dpi=200, label="R_in", xlabel="temperature in Kelvin", ylabel="Radiation in W/m²")
    plot!(temperatures, map(T -> R_out(T), temperatures), label="R_out")
end
gif(anim, "./MaKli/exercises/graphics/4a.gif", fps=30)

# unique intersection if R_in(290, Q(t)) < R_out(290)

# (b)
temperatures = (-30:40) .+ celsius_in_kelvin

function foo(temperatures)
    t_end1 = 3600 * 24 * 365 * 1 # 1 Jahr
    t_span1 = (0.0, t_end1)
    h1 = Int(ceil(t_end1 / 1e5))
    solution1 = zeros(eachindex(temperatures))

    for i in eachindex(temperatures)
        t, sol = explicit_euler_solve((T, t) -> u(first(T), t), [temperatures[i]], t_span1, h1)
        solution1[i] = sol[end]
    end

    t_end5 = 3600 * 24 * 365 * 5 # 5 Jahre
    t_span5 = (0.0, t_end5)
    h5 = Int(ceil(t_end5 / 1e5))
    solution5 = zeros(eachindex(temperatures))

    for i in eachindex(temperatures)
        t, sol = explicit_euler_solve((T, t) -> u(first(T), t), [temperatures[i]], t_span5, h5)
        solution5[i] = sol[end]
    end

    df = DataFrame("temperature in °C" => temperatures .- celsius_in_kelvin, "1 year" => solution1, "5 years" => solution5)
    df_long = stack(df, ["1 year", "5 years"], "temperature in °C", variable_name=:horizon, value_name=:temperature)
    return df_long
end

df = foo(temperatures)
Gadfly.plot(
    df,
    x="temperature in °C",
    y="temperature",
    color="horizon",
    Gadfly.Geom.point
)


# (c)
const T_star = (S₀ * (1 - α(T)))^(1 / 4)
g(T) = b * (1 - 0.4tanh(x / T_star)^6)
R_out(T) = g(T) * T^4

temperatures = (-30:0.1:30) .+ celsius_in_kelvin

anim = Plots.@animate for t in [(1 / ω_S) * (i / 100) for i in 1:100]
    plot(temperatures, map(T -> R_in(T, Q(t)), temperatures), dpi=200)
    plot!(temperatures, map(T -> R_out(T), temperatures))
end
gif(anim, "./MaKli/exercises/graphics/4c.gif", fps=30)

df = foo(temperatures)
Gadfly.plot(
    df,
    x="temperature in °C",
    y="temperature",
    color="horizon",
    Gadfly.Geom.point
)