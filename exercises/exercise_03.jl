using Pkg; Pkg.activate("./MaKli")
using MaKli
using DifferentialEquations
using Gadlfy

const α₀ = ...
const S₀ = ...
α(t) = α₀ * (1 - sin(ω * t) / 10)
S(t) = S₀ * (1 - 0.06 * sin(ω_hat * t))
const ω = 2π/(3600*24)
# const H = ...
# const ρ = ...
# const C = ...

f(T, t; H, ρ, C) = (S(t) * (1 - α(t)) / 4 - T^4 * (ε * σ)) / (H * ρ * C)


f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

