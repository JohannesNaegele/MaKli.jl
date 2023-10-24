using MaKli
using DataFrames
using CairoMakie

include("rahmstorf.jl")

function test_specs(params, param_seq, temperatures, salt, deviation=0.05) # use fixed salt concentration
    my_plots = []
    values::Vector{Float64} = [v for (k, v) in params]
    for i in eachindex(params)
        changed_values_upper = values .* [(i == j) ? 1.0 : 1.05 for j in eachindex(params)]
        changed_values_lower = values .* [(i == j) ? 1.0 : 1.0/1.05 for j in eachindex(params)]
        params_upper = Dict(param_seq[i] => changed_values_upper[i] for i in eachindex(param_seq))

        function foo(temperatures, salt)
            opt = (x, t) -> u(
                x;
                C=params_upper[:C],
                α=params_upper[:α],
                β=params_upper[:β],
                V1=params_upper[:V1],
                V2=params_upper[:V2],
                V3=params_upper[:V3],
                V4=params_upper[:V4],
                S0=params_upper[:S0],
                F1=params_upper[:F1],
                F2=params_upper[:F2],
                λ1=params_upper[:λ1],
                λ2=params_upper[:λ2],
                λ3=params_upper[:λ3],
                T1_star=params_upper[:T1_star],
                T2_star=params_upper[:T2_star],
                T3_star=params_upper[:T3_star]
            )

            t_end1 = 3600 * 24 * 365 * 1 # 1 Jahr
            t_span1 = (0.0, t_end1)
            h1 = Int(ceil(t_end1 / 1e5))
            solution1 = zeros(8, eachindex(temperatures))
        
            for i in eachindex(temperatures)
                t, sol = fast_explicit_euler_solve(
                    opt,
                    vcat([temperatures[i] for j in 1:4], [salt for j in 1:4]),
                    t_span1,
                    h1
                )
                solution1[:, i] = sol
            end
        
            t_end5 = 3600 * 24 * 365 * 5 # 5 Jahre
            t_span5 = (0.0, t_end5)
            h5 = Int(ceil(t_end5 / 1e5))
            solution5 = zeros(8, eachindex(temperatures))
        
            for i in eachindex(temperatures)
                t, sol = fast_explicit_euler_solve(
                    opt,
                    vcat([temperatures[i] for j in 1:4], [salt for j in 1:4]),
                    t_span5,
                    h5
                )
                solution5[:, i] = sol
            end
        
            return solution1, solution5
        end

        results = foo(temperatures, salt)
        println([temperatures[i][1] for i in eachindex(temperatures)])
        println([results[1][1, i] for i in eachindex(temperatures)])
        f = Figure(resolution = (600, 600));
        ax1 = Axis(f[1, 1], yticklabelcolor = :blue, ylabel = "Temperatur in Kelvin", xlabel="Anfangstemperatur", xreversed=true)
        ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, xreversed=true, ylabel="Salzkonzentration")
        hidespines!(ax2)
        hidexdecorations!(ax2)

        # watch out: if you have two axes you need the same x-axis length!
        lines!(ax1, [temperatures[i][1] for i in eachindex(temperatures)], [results[1][1, i] for i in eachindex(temperatures)], color = :blue)
        lines!(ax2, [temperatures[i][1] for i in eachindex(temperatures)], [results[2][4, i] for i in eachindex(temperatures)], color = :red)

        push!(my_plots, f)
    end
    return my_plots
end

# test_specs(params, (1:15) .+ celsius_in_kelvin, 33)
param_seq = [:C, :α, :β, :V1, :V2, :V3, :V4, :S0, :F1, :F2, :λ1, :λ2, :λ3, :T1_star, :T2_star, :T3_star]
@profview a = test_specs(parameters, param_seq, (1:15) .+ celsius_in_kelvin, 33)

(9.5 + 14.5 + 11.5 + 7 + 4.5) / (16 + 16 + 16 + 16 + 8)
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