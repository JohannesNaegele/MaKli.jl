using Pkg; Pkg.activate("./MaKli")
using DataFrames
using Gadfly

DataFrame(
    "Process" => [
        "collision of droplets during cloud formation",
        "formation of convection cells",
        "development of large-scale weather systems",
        "persistence of pressure distributions",
        "Southern Oscillation",
        "troposphere-stratosphere exchange"
    ],
    "Characteristic Time Scale" => [
        [1e-6, 1e-3],
        [1e4, 1e5],
        [1e4, 1e5],
        [1e6, 1e6],
        [1e7, 1e7],
        [1e7, 1e8],
    ],
    "Characteristic Spatial Scale" => [
        [1e-6, 1e-6],
        [1e2, 1e4],
        [1e6, 1e7],
        [1e6, 1e7],
        [1e7, 1e7],
        [0.0, Inf]
    ]
)