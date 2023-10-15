module MaKli

include("helpers/euler.jl")
include("helpers/equations.jl")

export explicit_euler_solve, implicit_euler_solve, @equations

end
