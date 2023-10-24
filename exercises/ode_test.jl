using DifferentialEquations

t_end = 20
prob = ODEProblem(
    (x, p, t) -> [0.01t^2 + x[2]^3, p[1] <= x[1] < p[1] + 0.5 ? 5.0 : 0.0],
    [1, 0],
    (0,t_end),
    [5.437]
)
sol = solve(prob, Rosenbrock23())
sol = solve(prob, TRBDF2())


h = 0.01
df = DataFrame(
    :time => 1:h:t_end,
    :x_1 => [sol(t)[1] for t in 1:h:t_end],
    :x_2 => [sol(t)[2] for t in 1:h:t_end]
)
df_long = stack(df, [:x_1, :x_2], :time, variable_name=:variable, value_name=:x)
plot(
    df_long,
    x=:time,
    y=:x,
    color=:variable,
    Geom.line
)