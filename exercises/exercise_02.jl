using DataFrames
using CairoMakie
using Pipe
using CSV
using Statistics
using GLM
using Optim
using NLsolve
using Impute
CairoMakie.activate!(type = "svg")

## Aufgabe 1

co2 = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_co2nat-noaa.csv", comment="#") |> DataFrame
temp = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_deutnat-noaa.csv", comment="#") |> DataFrame

df = outerjoin(co2, temp, on = :gas_ageBP => :ice_ageBP)
sort!(df, :gas_ageBP, rev=true)

# (a)

f = Figure(resolution = (1300, 600));
ax1 = Axis(f[1, 1], yticklabelcolor = :blue, ylabel = "CO2 in PPM", xlabel="Jahre vor 1999", xreversed=true, xtickformat = values -> ["$(Int(value)) Tsd." for value in values])
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, xreversed=true, ylabel="Temperaturdifferenz in Grad Celsius")
hidespines!(ax2)
hidexdecorations!(ax2)

# watch out: if you have two axes you need the same x-axis length!
lines!(ax1, df.gas_ageBP./1000, Impute.interp(df.CO2), color = :blue)
lines!(ax2, df.gas_ageBP./1000, Impute.interp(df.deltaTS), color = :red)

f
save("./exercises/co2_temp.svg", f)

# (b)
df_rounded = @pipe df |>
    transform!(_, :gas_ageBP => ByRow(x -> Int(round(x/1000) * 1000)) => :rounded_time) |>
    groupby(_, :rounded_time) |>
    combine(_,
        # :CO2 => mean ∘ skipmissing => :CO2,
        :deltaTS => mean ∘ skipmissing => :deltaTS_mean
    )
df_regr = @pipe df |>
    leftjoin(_, df_rounded, on = :rounded_time) |>
    _[!, [:gas_ageBP, :CO2, :deltaTS_mean]] |>
    dropmissing(_)
    # filter(r -> all(x -> !isnan(x), r), _)

cor(df_regr.CO2, df_regr.deltaTS_mean)
cor(exp.(df_regr.CO2), df_regr.deltaTS_mean)
cor(log.(df_regr.CO2), df_regr.deltaTS_mean)

# (c)
fm1 = @formula(deltaTS_mean ~ CO2)
fm2 = @formula(deltaTS_mean ~ exp(CO2))
# fm2 = @formula(log(deltaTS_mean+100) ~ CO2)

regr1 = lm(fm1, df_regr)
regr2 = lm(fm2, df_regr, dropcollinear=false)
r2(regr1)
r2(regr2)

predict(regr1, DataFrame(:CO2 => [233, 300, 420, 500]))
predict(regr2, DataFrame(:CO2 => [233, 300, 420, 500]))

# sanity check 1
f = Figure(resolution = (1300, 600));
ax1 = Axis(f[1, 1], yticklabelcolor = :blue, ylabel = "CO2 in PPM", xlabel="Temperaturdifferenz in Grad Celsius gerundet auf 1000 Jahre", xreversed=true, xtickformat = values -> ["$(Int(value)) Tsd." for value in values])

lines!(ax1, temp.ice_ageBP./1000, temp.deltaTS, color = :blue)
lines!(ax1, df_regr.gas_ageBP./1000, df_regr.deltaTS_mean, color = :red)
f

# sanity check 2
sort!(df_regr, :CO2)
f = Figure(resolution = (1300, 600));
ax1 = Axis(f[1, 1], yticklabelcolor = :blue, xlabel = "CO2 in PPM", ylabel = "Temperaturdifferenz in Grad Celsius gerundet auf 1000 Jahre")

scatter!(ax1, df_regr.CO2, df_regr.deltaTS_mean, color = :blue)
f
save("./exercises/co2_temp_cor.svg", f)

## Aufgabe 2
g(x) = 3 - (x * exp(x))/(exp(x) - 1)
sol = nlsolve(x -> g(first(x)), [1.0])
sol.zero

# ------------------------
# | experimenteller Teil |
# ------------------------

const h = 6.62607015e-34
const k_B = 1.380649e-23 # Boltzmann constant
const c = 299792458

function u(ν, T)
    first = (8π*h*ν^3)/(c^3)
    second = 1/(exp((h*ν) / (k_B*T)) - 1)
    return first * second
end

function u_maximizer(T)
    ν_initial = 2.821439 * k_B * T / h  # Wien's Law
    sol = optimize(ν -> -u(ν, T), ν_initial/2, ν_initial*2)  # Provide bounds around initial value
    return Optim.minimizer(sol)
end

# Convert frequencies to wavelengths and express in micrometers (µm)
wavelengths_µm = (c ./ u_maximizer.([1000, 1500, 2000])) * 1e6

# correct approach
println(wavelengths_µm)
c/(2.821439*k_B*5500/h)*1e6