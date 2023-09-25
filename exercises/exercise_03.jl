using Pkg; Pkg.activate("./MaKli")
using MaKli
using DifferentialEquations
using Gadfly
using DataFrames

const σ = 5.670374419e-8 # Stefan-Boltzmann constant
const ω = 2π/(3600*24)
const ω_hat = ω * 365
const α₀ = 0.3 # Albedo
const S₀ = 1361 # Solar constant
const ε = 0.62 # Emissionsgrad

H = 8007 # Skalenhöhe
ρ = 1.2 # Dichte
C = 1000

# (a)
α(t) = α₀ * (1 - sin(ω * t) / 10)
S(t) = S₀ * (1 - 0.06 * sin(ω_hat * t))

# t_end = (365*7*24*3600)
t_end = 40000000
f(T, t; H, ρ, C) = (S(t) * (1 - α(t)) / 4 - T^4 * (ε * σ)) / (H * ρ * C)
T_0 = 285.0

# idiomatic solution
t_span = (0.0, t_end)
prob = ODEProblem((T, p, t) -> f(T, t; H=H, ρ=ρ, C=C), T_0, t_span)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

# 
plot(
    x=1:100:t_end,
    y=[sol(t) for t in 1:100:t_end],
    Geom.line
)

# one week
plot(
    x=1:100:24*3600*7,
    y=[sol(t) for t in 1:100:24*3600*7],
    Geom.line
)

# (b)
H_W = 3
ρ_W = 1.5
C_W = 1000
const c = 0.6 # Wolkenbedeckung

celsius_in_kelvin = 273.15
T_0 = celsius_in_kelvin + 13
T_W_0 = celsius_in_kelvin - 45

# TODO: epsilon gleich?
ΔTΔt(T, T_W; H, ρ, C) = (S₀ * (1 - α₀) / 4 + c * σ * T_W^4 - T^4 * (ε * σ)) / (H * ρ * C)
ΔT_WΔt(T, T_W; H_W, ρ_W, C_W) = (c * ε * σ * T^4 - T_W^4 * σ * 2c) / (H_W * ρ_W * C_W)

# idiomatic solution
t_end = 3600*24*100 # 100 Tage
t_span = (0.0, t_end)
prob = ODEProblem((temp, p, t) -> [ΔTΔt(temp[1], temp[2]; H=H, ρ=ρ, C=C), ΔT_WΔt(temp[1], temp[2]; H_W=H_W, ρ_W=ρ_W, C_W=C_W)], [T_0, T_W_0], t_span)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

df = DataFrame(:t => 1:100:t_end, :T => [sol(t)[1] for t in 1:100:t_end], :T_W => [sol(t)[2] for t in 1:100:t_end])
df[!, :t] ./= 3600*24

df_long = stack(df, [:T, :T_W], :t, variable_name=:component, value_name=:temperature)
Gadfly.push_theme(:dark) # dark mode
plot(
    df_long,
    x=:t,
    y=:temperature,
    color=:component,
    Geom.line,
    Guide.xlabel("time in days")
)