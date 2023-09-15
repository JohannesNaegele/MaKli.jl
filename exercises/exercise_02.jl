using Pkg; Pkg.activate("./MaKli")
using DataFrames
using CairoMakie
using Pipe
using CSV
using Statistics
using GLM

co2 = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_co2nat-noaa.csv", comment="#") |> DataFrame
temp = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_deutnat-noaa.csv", comment="#") |> DataFrame

df = leftjoin(sort(co2, :gas_ageBP), sort(temp, :ice_ageBP), on = :gas_ageBP => :ice_ageBP)

# (a)

f = Figure();
ax1 = Axis(f[1, 1], yticklabelcolor = :blue, ylabel = "CO2 in PPM", xlabel="Years before 1999")
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, ylabel="temperature difference in Kelvin")
hidespines!(ax2)
hidexdecorations!(ax2)

lines!(ax1, co2.gas_ageBP, co2.CO2, color = :blue)
lines!(ax2, temp.ice_ageBP, temp.deltaTS, color = :red)


f

# 1999 - gas_ageBP

# (b)
df_regr = innerjoin(co2, temp, on = :gas_ageBP => :ice_ageBP)

cor(df_regr.co2, df_regr.temp)

# (c)
fm1 = @formula(temp ~ co2)
fm2 = @formula(temp ~ exp(co2))

regr1 = lm(df_regr, fm1)
regr2 = lm(df_regr, fm2)