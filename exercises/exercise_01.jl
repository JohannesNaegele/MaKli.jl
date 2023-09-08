using Pkg; Pkg.activate("./MaKli")
using DataFrames
using Gadfly
using Pipe
import Cairo, Fontconfig

Gadfly.push_theme(:default) # dark mode
# Gadfly.push_theme(style(background_color=nothing))

include("components.jl")

dfs = [atmosphere, hydrosphere, cryosphere, land_surface, biosphere]
ids = ["atmosphere", "hydrosphere", "cryosphere", "land surface", "biosphere"]

for i in eachindex(dfs)
    transform!(dfs[i], :Process => (x -> ids[i]) => :component)
end

df = vcat(atmosphere, hydrosphere, cryosphere, land_surface, biosphere)

margin(x) = x[1] == x[2] ? [log10(x[1]) - 0.5, log10(x[1]) + 0.5].^10 : x
function bound(x)
    highest = maximum(filter(!=(Inf), x))
    lowest = minimum(filter(!=(0.0), x))
    for i in eachindex(x)
        if x[i] == Inf
            x[i] = 1e5^ceil(log(1e5, highest))
        elseif x[i] == 0.0
            x[i] = 1e5^floor(log(1e5, lowest))
        end
    end
    return x
end

# @pipe df |>
#     transform!(
#         _,
#         ["Characteristic Time Scale", "Characteristic Spatial Scale"] => ByRow((x, y) -> [margin(x), margin(y)]) => ["Characteristic Time Scale", "Characteristic Spatial Scale"]
#     ) |>
#     transform!(
#         _,
#         ["Characteristic Time Scale", "Characteristic Spatial Scale"] => ByRow((x, y) -> [x..., y...]) => [:x_min, :x_max, :y_min, :y_max]
#     ) |>
#     transform!(
#         _,
#         [:x_min, :x_max, :y_min, :y_max] => (u, x, y, z) -> bound.([u, x, y, z]) => [:x_min, :x_max, :y_min, :y_max]
#     )

# plot(
#     atmosphere,
#     x="Characteristic Time Scale",
#     y="Characteristic Spatial Scale",
#     # color=:id, 
#     alpha=[0.3],
#     linestyle=[:dash],
#     Geom.polygon(fill=true),
#     Scale.color_discrete,
#    Theme(line_width=2pt,
#    lowlight_color=identity,
#    discrete_highlight_color=identity)
# )

plot(
    df,
    xmin=:x_min,
    xmax=:x_max,
    ymin=:y_min,
    ymax=:y_max,
    color=:component,
    alpha=[0.8],
    Geom.rect,
    Scale.x_log10, Scale.y_log10,
    Guide.xlabel("Characteristic Time Scale"),
    Guide.ylabel("Characteristic Spatial Scale")
)

p = plot(
    df,
    xmin=:x_min,
    xmax=:x_max,
    ymin=:y_min,
    ymax=:y_max,
    color=:Process,
    alpha=[0.5],
    ygroup=:component,
    Geom.subplot_grid(Geom.rect),
    Scale.x_log10, Scale.y_log10,
    Guide.xlabel("Characteristic Time Scale in s"),
    Guide.ylabel("Characteristic Spatial Scale in m")
)

draw(PDF("components.pdf", 20cm, 20cm), p)

# TODO: fix ^ Unicode

# @pipe df |> subset(_, :component => ByRow(==("atmosphere"))) |>
#     plot(
#         _,
#         xmin=:x_min,
#         xmax=:x_max,
#         ymin=:y_min,
#         ymax=:y_max,
#         color=:Process,
#         alpha=[0.5],
#         Geom.rect,
#         Scale.x_log10, Scale.y_log10,
#         Guide.xlabel("Characteristic Time Scale"),
#         Guide.ylabel("Characteristic Spatial Scale")
#     )