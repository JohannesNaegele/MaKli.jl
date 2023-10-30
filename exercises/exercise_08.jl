using MaKli
using Plots
using Random

const T = 100
const n = Int(1e6)
function u!(dx, x, dw)
    dx[1] = atan(1 / (1 + x[1]^2)) + sin(x[1]) * dw
end
# [0; sqrt(T / n) .* cumsum(randn(n))]
Random.seed!(1)

plotly()
p = plot();
for i in 1:4
    t_eul, sol_eul = explicit_euler_solve!(
        (du, u, t) -> (u!(du, u, randn())),
        [1.0],
        (0, T),
        T/n
    )
    plot!(collect(t_eul), vec(sol_eul))
end
p
savefig(p)