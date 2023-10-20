using NLsolve

function explicit_euler_solve(f, u0, t_span, h)
    t0, t_end = t_span
    t = t0:h:t_end
    u = zeros(length(u0), length(t))
    u[:, 1] = u0
    
    for i = 1:length(t)-1
        u[:, i+1] .= u[:, i] .+ h * f(u[:, i], t[i])
    end
    
    return t, u
end

function explicit_euler_solve!(f!, u0, t_span, h)
    t0, t_end = t_span
    t = t0:h:t_end
    u = zeros(length(u0), length(t))
    du = zeros(length(u0))
    u[:, 1] = u0
    
    for i = 1:length(t)-1
        f!(du, u[:, i], t[i])
        u[:, i+1] .= u[:, i] .+ h * du
    end
    
    return t, u
end

function fast_explicit_euler_solve(f, u0, t_span, h)
    t0, t_end = t_span
    t = t0:h:t_end
    u = deepcopy(u0)
    
    for i = 1:length(t)-1
        u = u .+ h * f(u, t[i])
    end
    
    return t, u
end

function implicit_euler_solve(f, u0, t_span, h)
    t0, t_end = t_span
    t = t0:h:t_end
    u = zeros(length(u0), length(t))
    u[:, 1] = u0
    
    for i = 1:length(t)-1    
        function g!(u_next)
            u_next .- u[:, i] .- h * f(u_next, t[i+1])
        end

        result = nlsolve(g!, u[:, i])
        u[:, i+1] = result.zero
    end
    
    return t, u
end
