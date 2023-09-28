using Pkg; Pkg.activate("./MaKli")
using MaKli
using DifferentialEquations
using Gadfly
using DataFrames
import Cairo, Fontconfig

const σ = 5.670374419e-8 # Stefan-Boltzmann constant
const ω = 2π / (3600 * 24)
const ω_hat = ω / 365
const α₀ = 0.3 # Albedo
const S₀ = 1361 # Solar constant
const ε = 0.62 # Emissionsgrad

H = 8007 # Skalenhöhe
ρ = 1.2 # Dichte
C = 1000

# (a)
α(t) = α₀ * (1 - sin(ω * t) / 10)
S(t) = S₀ * (1 - 0.06 * sin(ω_hat * t))

day_in_seconds = 24 * 3600
t_end = (1000 * day_in_seconds) # 1000 days
# t_end = 40000000
f(T, t; H, ρ, C) = (S(t) * (1 - α(t)) / 4 - T^4 * (ε * σ)) / (H * ρ * C)
T_0 = 285.0
t_span = (0.0, t_end)
h = 60 # 1 min

# idiomatic solution
prob = ODEProblem((T, p, t) -> f(T, t; H=H, ρ=ρ, C=C), T_0, t_span)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

plot(
    x=1:h:t_end,
    y=[sol(t) for t in 1:h:t_end],
    Geom.line
)

# custom solution
t, u = explicit_euler_solve((T, t) -> f(first(T), t; H=H, ρ=ρ, C=C), [T_0], t_span, h)

p1 = plot(
    x=t ./ day_in_seconds,
    y=u,
    Geom.line,
    Guide.title("Solution with explicit Euler"),
    Guide.xlabel("time in days"),
    Guide.ylabel("temperature in Kelvin")
)

# one week
p2 = plot(
    x=(t./day_in_seconds)[begin:(60*24*8)],
    y=u[begin:(60*24*8)],
    Geom.line,
    Guide.title("Solution with explicit Euler"),
    Guide.xlabel("time in days"),
    Guide.ylabel("temperature in Kelvin")
)

t, u = implicit_euler_solve((T, t) -> f(first(T), t; H=H, ρ=ρ, C=C), [T_0], t_span, h)

p3 = plot(
    x=t ./ day_in_seconds,
    y=u,
    Geom.line,
    Guide.title("Solution with implicit Euler"),
    Guide.xlabel("time in days"),
    Guide.ylabel("temperature in Kelvin")
)

# one week
p4 = plot(
    x=(t./day_in_seconds)[begin:(60*24*8)],
    y=u[begin:(60*24*8)],
    Geom.line,
    Guide.title("Solution with implicit Euler"),
    Guide.xlabel("time in days"),
    Guide.ylabel("temperature in Kelvin")
)

# (b)
H_W = 3
ρ_W = 1.5
C_W = 1000
const c = 0.6 # Wolkenbedeckung

celsius_in_kelvin = 273.15
T_0 = celsius_in_kelvin + 13
T_W_0 = celsius_in_kelvin - 45
t_end = 3600 * 24 * 100 # 100 Tage
t_span = (0.0, t_end)
h = Int(ceil(t_end / 1e5))

# TODO: epsilon gleich?
ΔTΔt(T, T_W; H, ρ, C) = (S₀ * (1 - α₀) / 4 + c * σ * T_W^4 - (ε * σ) * T^4) / (H * ρ * C)
ΔT_WΔt(T, T_W; H_W, ρ_W, C_W) = (c * ε * σ * T^4 -  2c * σ * T_W^4) / (H_W * ρ_W * C_W)

# idiomatic solution
prob = ODEProblem((temp, p, t) -> [ΔTΔt(temp[1], temp[2]; H=H, ρ=ρ, C=C), ΔT_WΔt(temp[1], temp[2]; H_W=H_W, ρ_W=ρ_W, C_W=C_W)], [T_0, T_W_0], t_span)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

df = DataFrame(:t => 1:100:t_end, :T => [sol(t)[1] for t in 1:100:t_end], :T_W => [sol(t)[2] for t in 1:100:t_end])
df[!, :t] ./= 3600 * 24

df_long = stack(df, [:T, :T_W], :t, variable_name=:component, value_name=:temperature)

plot(
    df_long,
    x=:t,
    y=:temperature,
    color=:component,
    Geom.line,
    Guide.xlabel("time in days")
)

# explicit euler
t, u = explicit_euler_solve((temp, t) -> [ΔTΔt(temp[1], temp[2]; H=H, ρ=ρ, C=C), ΔT_WΔt(temp[1], temp[2]; H_W=H_W, ρ_W=ρ_W, C_W=C_W)], [T_0, T_W_0], t_span, h)
df = DataFrame(:t => t, :T => u[1, :], :T_W => u[2, :])
df[!, :t] ./= 3600 * 24

df_long = stack(df, [:T, :T_W], :t, variable_name=:component, value_name=:temperature)

p5 = plot(
    df_long,
    x=:t,
    y=:temperature,
    color=:component,
    Geom.line,
    Guide.xlabel("time in days"),
    Guide.ylabel("temperature in Kelvin"),
    Guide.title("Solution with explicit Euler")
)

# implicit euler
t, u = implicit_euler_solve((temp, t) -> [ΔTΔt(temp[1], temp[2]; H=H, ρ=ρ, C=C), ΔT_WΔt(temp[1], temp[2]; H_W=H_W, ρ_W=ρ_W, C_W=C_W)], [T_0, T_W_0], t_span, h)
df = DataFrame(:t => t, :T => u[1, :], :T_W => u[2, :])
df[!, :t] ./= 3600 * 24

df_long = stack(df, [:T, :T_W], :t, variable_name=:component, value_name=:temperature)

p6 = plot(
    df_long,
    x=:t,
    y=:temperature,
    color=:component,
    Geom.line,
    Guide.xlabel("time in days"),
    Guide.ylabel("temperature in Kelvin"),
    Guide.title("Solution with implicit Euler")
)

draw(PDF("./MaKli/exercises/graphics/3a_euler_explicit_1.pdf", 14cm, 14cm), p1)
draw(PDF("./MaKli/exercises/graphics/3a_euler_explicit_2.pdf", 14cm, 14cm), p2)
draw(PDF("./MaKli/exercises/graphics/3a_euler_implicit_1.pdf", 14cm, 14cm), p3)
draw(PDF("./MaKli/exercises/graphics/3a_euler_implicit_2.pdf", 14cm, 14cm), p4)
draw(PDF("./MaKli/exercises/graphics/3b_euler_explicit.pdf", 14cm, 14cm), p5)
draw(PDF("./MaKli/exercises/graphics/3b_euler_implicit.pdf", 14cm, 14cm), p6)