using MaKli
using Gadfly
using DataFrames

params = @equations begin
    λ1 = 3000
    λ2 = 3000
    λ3 = 1000
    S0 = 35
    T1_star = 279.6
    T2_star = 275.7
    T3_star = 284.7
    V1 = 1.1
    V2 = 0.4
    V3 = 0.68
    V4 = 0.05
end

function Δtemp(T)
    return [
        (M(t) / V1) * (T[4] - T[1]) + λ1 * (T1_star - T[1]),
        (M(t) / V2) * (T[3] - T[2]) + λ2 * (T2_star - T[2]),
        (M(t) / V3) * (T[1] - T[3]) + λ3 * (T3_star - T[3]),
        (M(t) / V4) * (T[2] - T[4])
    ]
end

function Δsalt(S)
    return [
        (M(t) / V1) * (S[4] - S[1]) + S0 * F1 / V1
        (M(t) / V2) * (S[3] - S[2]) + S0 * F2 / V2
        (M(t) / V3) * (S[1] - S[3]) + S0 * (F2 - F1) / V3
        (M(t) / V4) * (S[2] - S[4])
    ]
end

function u(x)
    return vcat(Δtemp(x[1:4]), Δsalt([1:4]))
end


## ODE test

# using DifferentialEquations

# t_end = 20
# prob = ODEProblem(
#     (x, p, t) -> [0.01t^2 + x[2]^3, p[1] <= x[1] < p[1] + 0.5 ? 5.0 : 0.0],
#     [1, 0],
#     (0,t_end),
#     [5.437]
# )
# sol = solve(prob, Rosenbrock23())
# sol = solve(prob, TRBDF2())


# h = 0.01
# df = DataFrame(
#     :time => 1:h:t_end,
#     :x_1 => [sol(t)[1] for t in 1:h:t_end],
#     :x_2 => [sol(t)[2] for t in 1:h:t_end]
# )
# df_long = stack(df, [:x_1, :x_2], :time, variable_name=:variable, value_name=:x)
# plot(
#     df_long,
#     x=:time,
#     y=:x,
#     color=:variable,
#     Geom.line
# )