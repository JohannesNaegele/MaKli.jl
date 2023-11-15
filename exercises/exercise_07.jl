using MaKli
using DataFrames
using Plots
using Distances
using DifferentialEquations

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
h = t_end / n_steps

# Parameter choice
x0 = [10, 13, 34]

# Euler solution
t, sol = explicit_euler_solve!(
    (du, u, t) -> (lorenz!(du, u, σ=9, R=23, B=3)),
    x0,
    t_span,
    h
)

plot(
    t,
    [sol[1, :], sol[2, :], sol[3, :]],
    layout=(3, 1),
    legend=false,
    title="Explicit Euler"
)
plot(
    sol[1, :], sol[2, :], sol[3, :],
    title="Explicit Euler",
    label="x,y,z"
)
last_third = Int(floor(length(sol[1, :]) * 2 / 3)):length(sol[1, :])
plot(
    sol[1, last_third], sol[2, last_third], sol[3, last_third],
    title="Explicit Euler, last third of iterations",
    label="x,y,z"
)

# Runge-Kutta solution
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    x0,
    t_span,
    DP5(),
)
sol = solve(prob)
last_third = Int(floor(length(sol.t) * 2 / 3)):length(sol.t)
plot(sol.t, [sol[1, :], sol[2, :], sol[3, :]], layout=(3, 1), legend=false)
plot(sol, idxs=(1, 2, 3), title="Dormand-Prince's explicit 5/4 Runge-Kutta method")
plot(
    sol[1, last_third], sol[2, last_third], sol[3, last_third],
    title="Dormand-Prince's explicit 5/4 Runge-Kutta method",
    label="x,y,z"
)

# With reltol
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    x0,
    t_span,
    DP5(),
    reltol=1e-7
)
sol = solve(prob)
plot(sol.t, [sol[1, :], sol[2, :], sol[3, :]], layout=(3, 1), legend=false)
plot(sol, idxs=(1, 2, 3), title="Dormand-Prince's explicit 5/4 Runge-Kutta method\nwith reltol")
last_third = Int(floor(length(sol.t) * 2 / 3)):length(sol.t)
plot(
    sol[1, last_third], sol[2, last_third], sol[3, last_third],
    title="Dormand-Prince's explicit 5/4 Runge-Kutta method\nwith reltol",
    label="x,y,z"
)

# (d)
x0 = [1.0, 1.0, 1.0]
test_grid = 10:1:110
function test_fixed(x0, grid, last=20, solver=DP5(), rt=false)
    p = plot()
    for r in grid
        prob = ODEProblem(
            (du, u, p, t) -> lorenz!(du, u, σ=10, R=r, B=8/3),
            x0,
            t_span,
            solver,
            reltol=rt
        )
        sol = solve(prob)
        plot!([r for i in 1:last], sol[1, (end-last+1):end], seriestype=:scatter, legend=false)
    end
    return p
end

test_fixed(x0, test_grid);
title!("Letzte 20 Iterationen für verschiedene R\nmit 5/4 Runge-Kutta und reltol=1e-10");
savefig("./exercises/graphics/7d_1.pdf")
test_fixed(x0, 20:0.1:25, 20, DP5(), 1e-10);
title!("Letzte 20 Iterationen für verschiedene R\nmit 5/4 Runge-Kutta und reltol=1e-10");
savefig(savefig("./exercises/graphics/7d_2.pdf"))

function periodic(sol, tol)
    for i in eachindex(sol)[1:end-1]
        for j in eachindex(sol)[(i+1):end]
            if euclidean(sol[i], sol[j]) < tol
                return true
            end
        end
    end
    return false
end

function test_periodic(x0, grid, solver=DP5(), rt=false)
    p = plot()
    period = zeros(length(grid))
    for i in eachindex(grid)
        prob = ODEProblem(
            (du, u, p, t) -> lorenz!(du, u, σ=10, R=grid[i], B=8/3),
            x0,
            t_span,
            solver,
            reltol=rt
        )
        sol = solve(prob)
        if periodic(sol, 1e-4)
            period[i] = 1.0
        else
            period[i] = 0.0
        end
    end
    plot!(collect(grid), period, legend=false)
    return p
end

test_periodic(x0, 10:1:110, DP5(), 1e-10);
title!("Test auf Periodizität mit tol=1e-4\nund 5/4 Runge-Kutta und reltol=1e-10");
savefig(savefig("./exercises/graphics/7d_3.pdf"))
test_periodic(x0, 99:0.05:101, DP5(), 1e-10);
title!("Test auf Periodizität mit tol=1e-4\nund 5/4 Runge-Kutta und reltol=1e-10");
savefig(savefig("./exercises/graphics/7d_4.pdf"))

# (10, 120, 0.2)
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=10, R=120, B=0.2),
    x0,
    t_span,
    DP5(),
    reltol=1e-10
)
sol = solve(prob)
plot(sol.t, [sol[1, :], sol[2, :], sol[3, :]], layout=(3, 1), legend=false,
    plot_title="5/4 Runge-Kutta method\nwith reltol=1e-10"
);
savefig(savefig("./exercises/graphics/7d_5.pdf"))
plot(sol, idxs=(1, 2, 3), title="5/4 Runge-Kutta method\nwith reltol=1e-10", label="x,y,z");
last_third = Int(floor(length(sol.t) * 2 / 3)):length(sol.t)
savefig(savefig("./exercises/graphics/7d_6.pdf"))
plot(
    sol[1, last_third], sol[2, last_third], sol[3, last_third],
    title="5/4 Runge-Kutta method\nwith reltol=1e-10, last third",
    label="x,y,z"
);
savefig(savefig("./exercises/graphics/7d_7.pdf"))

# (10, 118, 0.1)
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=10, R=118, B=0.1),
    x0,
    t_span,
    DP5(),
    reltol=1e-10
)
sol = solve(prob)
plot(sol.t, [sol[1, :], sol[2, :], sol[3, :]], layout=(3, 1), legend=false,
    plot_title="5/4 Runge-Kutta method\nwith reltol=1e-10"
);
savefig(savefig("./exercises/graphics/7d_8.pdf"))
plot(sol, idxs=(1, 2, 3), title="5/4 Runge-Kutta method\nwith reltol=1e-10", label="x,y,z");
savefig(savefig("./exercises/graphics/7d_9.pdf"))
last_third = Int(floor(length(sol.t) * 2 / 3)):length(sol.t)
plot(
    sol[1, last_third], sol[2, last_third], sol[3, last_third],
    title="5/4 Runge-Kutta method\nwith reltol=1e-10, last third",
    label="x,y,z"
);
savefig(savefig("./exercises/graphics/7d_10.pdf"))