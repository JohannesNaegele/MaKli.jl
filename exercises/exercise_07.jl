using MaKli
using DataFrames
using Plots
import CairoMakie
using DifferentialEquations

const celsius_in_kelvin = 273.15

# (c)
# Lorenz attractor
function lorenz!(dx, x; σ, R, B)
    dx[1] = -σ * x[1] + σ * x[2]
    dx[2] = -x[1] * x[3] + R * x[1] - x[2]
    dx[3] = x[1] * x[2] - B * x[3]
end

# Solver options
n_steps = 1e6
t_end = 200
t_span = (0.0, t_end)
h1 = Int(ceil(t_end / n_steps))

# Parameter choice
# σ, R, B = 9, 23, 3
x0 = [10, 13, 34]

# Euler solution
t, sol = explicit_euler_solve(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    vcat([temperatures[i] for j in 1:4], [salt for j in 1:4]),
    t_span1,
    h1
)
plot(t, sol, title="")

# Runge-Kutta solution
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    x0,
    t_span,
    DP5(),
)
sol = solve(prob)
plot(sol, idxs = (1, 2, 3), title="Dormand-Prince's explicit 5/4 Runge-Kutta method")

# With reltol
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    x0,
    t_span,
    DP5(),
    reltol=1e-7
)
sol = solve(prob)
plot(sol, idxs = (1, 2, 3), title="Dormand-Prince's explicit 5/4 Runge-Kutta method\nwith reltol")
# (d)

# (10, 20, 8/3)
# (10, 100, 8/3)
# (10, 28, 8/3)

# (10, 120, 0.2)
# (10, 118, 0.1)

# Tsit5()