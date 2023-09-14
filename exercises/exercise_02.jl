using Pkg; Pkg.activate("./MaKli")
using DataFrames
using Gadfly
using Pipe
using CSV

co2 = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_co2nat-noaa.csv", comment="#") |> DataFrame
temp = CSV.File("./data/ncei.noaa.gov_pub_data_paleo_icecore_antarctica_vostok_deutnat-noaa.csv", comment="#") |> DataFrame

