const s_year = 365 * 24 * 60 * 60
const celsius_in_kelvin = 273.15
const Γ = 7.3e8
const ρ0 = 1025
const c = 4000

parameters = @equations begin
    λ1 = Γ / (3000 * ρ0 * c)
    λ2 = Γ / (3000 * ρ0 * c)
    λ3 = Γ / (1000 * ρ0 * c)
    S0 = 35
    T1_star = 279.6
    T2_star = 275.7
    T3_star = 284.7
    V1 = 1.1
    V2 = 0.4
    V3 = 0.68
    V4 = 0.05
    C = 25.4
    α = 1.7e-4
    β = 8.0e-4
    F1 = (1.0e17 / s_year)
    F2 = (1.0e17 / s_year)
end

function M(T, S, C, α, β)
    return C * (-α + (T[2] - T[1]) + β * (S[2] - S[1]))
end

function Δtemp(T, S; C, α, β, V1, V2, V3, V4, λ1, λ2, λ3, T1_star, T2_star, T3_star)
    return [
        (M(T, S, C, α, β) / V1) * (T[4] - T[1]) + λ1 * (T1_star - T[1]),
        (M(T, S, C, α, β) / V2) * (T[3] - T[2]) + λ2 * (T2_star - T[2]),
        (M(T, S, C, α, β) / V3) * (T[1] - T[3]) + λ3 * (T3_star - T[3]),
        (M(T, S, C, α, β) / V4) * (T[2] - T[4])
    ]
end

function Δsalt(T, S; C, α, β, V1, V2, V3, V4, S0, F1, F2)
    return [
        (M(T, S, C, α, β) / V1) * (S[4] - S[1]) + S0 * F1 / V1,
        (M(T, S, C, α, β) / V2) * (S[3] - S[2]) + S0 * F2 / V2,
        (M(T, S, C, α, β) / V3) * (S[1] - S[3]) + S0 * (F2 - F1) / V3,
        (M(T, S, C, α, β) / V4) * (S[2] - S[4])
    ]
end

function u(x; C, α, β, V1, V2, V3, V4, S0, F1, F2, λ1, λ2, λ3, T1_star, T2_star, T3_star)
    return vcat(
        Δtemp(x[1:4], x[5:8], C=C, α=α, β=β, V1=V1, V2=V2, V3=V3, V4=V4, λ1=λ1, λ2=λ2, λ3=λ3, T1_star=T1_star, T2_star=T2_star, T3_star=T3_star),
        Δsalt(x[1:4], x[5:8], C=C, α=α, β=β, V1=V1, V2=V2, V3=V3, V4=V4, S0=S0, F1=F1, F2=F2)
    )
end