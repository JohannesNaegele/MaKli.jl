using MaKli
using DifferentialEquations
using Gadfly
using DataFrames
import Cairo, Fontconfig
using GLM
using CSV
using Pipe
Gadfly.push_theme(:dark)

# (a)
df = CSV.File("./data/webbook.nist.gov_cgi_fluid.csv", delim="\t") |> DataFrame

@pipe df |>
    subset!(_, "Temperature (C)" => ByRow(x -> 0.0 <= x <= 50.0))

df[!, "ρ"] = df[!, "Density (g/ml)"]
max_index = argmax(df[!, "ρ"])
df[!, "T"] = df[!, "Temperature (C)"] .- df[max_index, "Temperature (C)"]
formula = @formula(ρ ~ T + T^2 + T^3 + T^4)
regr = lm(formula, df)
r2(regr)
df[!, "Density fitted (g/ml)"] = predict(regr, df)

# relative error
maximum(abs.(df[!, "Density fitted (g/ml)"] .- df[!, "Density (g/ml)"]) ./ df[!, "Density (g/ml)"])

df_long = @pipe df |>
    stack(_, ["Density (g/ml)", "Density fitted (g/ml)"], "Temperature (C)", variable_name="Series type", value_name="g/ml")

p1 = plot(
    df_long,
    x="Temperature (C)",
    y="g/ml",
    color="Series type",
    Geom.line
);

df_long = @pipe df |>
    subset(_, "Temperature (C)" => ByRow(x -> 0.0 <= x <= 10.0)) |> 
    stack(_, ["Density (g/ml)", "Density fitted (g/ml)"], "Temperature (C)", variable_name="Series type", value_name="g/ml")

p2 = plot(
    df_long,
    x="Temperature (C)",
    y="g/ml",
    color="Series type",
    Geom.line
);

# (b)

draw(PDF("./exercises/graphics/5a_1.pdf", 14cm, 10cm), p1)
draw(PDF("./exercises/graphics/5a_2.pdf", 14cm, 10cm), p2)
draw(SVG("./exercises/graphics/5a_1.svg", 14cm, 10cm), p1)
draw(SVG("./exercises/graphics/5a_2.svg", 14cm, 10cm), p2)