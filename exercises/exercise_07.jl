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

# Initial values
x0 = [10, 13, 34]

# Euler solution
t_eul, sol_eul = explicit_euler_solve!(
    (du, u, t) -> (lorenz!(du, u, σ=9, R=23, B=3)),
    x0,
    t_span,
    h
)
last_third_eul = Int(floor(length(sol_eul[1, :]) * 2 / 3)):length(sol_eul[1, :])

# Runge-Kutta solution
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    x0,
    t_span,
    DP5(),
)
sol_45 = solve(prob)
last_third_45 = Int(floor(length(sol_45.t) * 2 / 3)):length(sol_45.t)

# With reltol
prob = ODEProblem(
    (du, u, p, t) -> lorenz!(du, u, σ=9, R=23, B=3),
    x0,
    t_span,
    DP5(),
    reltol=1e-7
)
sol_reltol = solve(prob)
last_third_reltol = Int(floor(length(sol_reltol.t) * 2 / 3)):length(sol_reltol.t)

plot(
    t_eul,
    [sol_eul[1, :], sol_eul[2, :], sol_eul[3, :]],
    layout=(3, 1),
    legend=false,
    plot_title="Explicit Euler"
)
plot(
    sol_eul[1, :], sol_eul[2, :], sol_eul[3, :],
    title="Explicit Euler",
    label="x,y,z"
)
plot(
    sol_eul[1, last_third_eul], sol_eul[2, last_third_eul], sol_eul[3, last_third_eul],
    title="Explicit Euler, last third of iterations",
    label="x,y,z"
)
plot!(sol_45.t, [sol_45[1, :], sol_45[2, :], sol_45[3, :]], layout=(3, 1), legend=false)
plot(sol_45, idxs=(1, 2, 3), title="Dormand-Prince's explicit 5/4 Runge-Kutta method")
plot(
    sol_45[1, last_third_45], sol_45[2, last_third_45], sol_45[3, last_third_45],
    title="Dormand-Prince's explicit 5/4 Runge-Kutta method",
    label="x,y,z"
)
plot!(sol_reltol.t, [sol_reltol[1, :], sol_reltol[2, :], sol_reltol[3, :]], layout=(3, 1), legend=false)
plot(sol_reltol, idxs=(1, 2, 3), title="Dormand-Prince's explicit 5/4 Runge-Kutta method\nwith reltol")
plot(
    sol_reltol[1, last_third_reltol], sol_reltol[2, last_third_reltol], sol_reltol[3, last_third_reltol],
    title="Dormand-Prince's explicit 5/4 Runge-Kutta method\nwith reltol=1e-7",
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