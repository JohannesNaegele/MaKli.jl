using Pkg; Pkg.activate("./MaKli")
using DataFrames
using Gadfly
using Pipe

include("components.jl")

dfs = [atmosphere, hydrosphere, cryosphere, land_surface, biosphere]
ids = ["atmosphere", "hydrosphere", "cryosphere", "land_surface", "biosphere"]

for i in eachindex(dfs)
    transform!(dfs[i], :Process => (x -> ids[i]) => :component)
end

df = vcat(atmosphere, hydrosphere, cryosphere, land_surface, biosphere)

# @pipe atmosphere |>
#     transform!(
#         _,
#         ["Characteristic Time Scale", "Characteristic Spatial Scale"] => ByRow()
#     )

process_repeated = repeat(atmosphere[!, "Process"], inner=4)

# Flatten vector columns
time_scale_flattened = vcat(atmosphere[!, "Characteristic Time Scale"]...)

plot(
    atmosphere,
    x="Characteristic Time Scale",
    y="Characteristic Spatial Scale",
    # color=:id, 
    alpha=[0.3],
    linestyle=[:dash],
    Geom.polygon(fill=true),
    Scale.color_discrete,
   Theme(line_width=2pt,
   lowlight_color=identity,
   discrete_highlight_color=identity)
)

reduce(vcat, [DataFrame(x=[1, 5, 5, 1].+d, y=[14, 14, 10, 10].-d, id=d+1) for d in 0:9])
[DataFrame(x=[1, 5, 5, 1].+d, y=[14, 14, 10, 10].-d, id=d+1) for d in 0:9]