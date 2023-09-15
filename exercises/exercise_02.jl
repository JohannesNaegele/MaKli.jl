using Pkg; Pkg.activate("./MaKli")
using DataFrames
using CairoMakie
using Pipe
using CSV
using Statistics
using GLM
CairoMakie.activate!(type = "svg")

co2 = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_co2nat-noaa.csv", comment="#") |> DataFrame
temp = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_deutnat-noaa.csv", comment="#") |> DataFrame

df = outerjoin(sort(co2, :gas_ageBP), sort(temp, :ice_ageBP), on = :gas_ageBP => :ice_ageBP)

# (a)

f = Figure(resolution = (1300, 600));
ax1 = Axis(f[1, 1], yticklabelcolor = :blue, ylabel = "CO2 in PPM", xlabel="Jahre vor 1999", xreversed=true, xtickformat = values -> ["$(Int(value)) Tsd." for value in values])
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, xreversed=true, ylabel="Temperaturdifferenz in Grad Celsius")
hidespines!(ax2)
hidexdecorations!(ax2)

lines!(ax1, co2.gas_ageBP./1000, co2.CO2, color = :blue)
lines!(ax2, temp.ice_ageBP./1000, temp.deltaTS, color = :red)

f
save("./exercises/co2_temp.svg", f)

# TODO: invert x-axis

# (b)
df_rounded = @pipe df |>
    transform(_, :gas_ageBP => ByRow(x -> round(x/1000) * 1000) => :gas_ageBP)
# round
df_rounded = outerjoin(roundco2, temp, on = :gas_ageBP => :ice_ageBP)
df_regr = @pipe df_rounded[!, [:gas_ageBP, :CO2, :deltaTS]] |> dropmissing

cor(df_regr.co2, df_regr.temp)

# (c)
fm1 = @formula(temp ~ co2)
fm2 = @formula(temp ~ exp(co2))

regr1 = lm(df_regr, fm1)
regr2 = lm(df_regr, fm2)