atmosphere = DataFrame(
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

hydrosphere = DataFrame(
    "Process" => [
        "gas exchange atmosphere-ocean",
        "deep water formation",
        "meso-scale oceanic gyres",
        "propagation of Rossby waves",
        "El Niño",
        "turnover of deep water"
    ],
    "Characteristic Time Scale" => [
        [1e-3, 1e6],
        [1e4, 1e6],
        [1e6, 1e7],
        [1e7, 1e7],
        [1e7, 1e8],
        [1e9, 1e10],
    ],
    "Characteristic Spatial Scale" => [
        [1e-6, 1e3],
        [1e4, 1e5],
        [1e4, 1e5],
        [1e7, 1e7],
        [1e7, 1e7],
        [0.0, Inf]
    ]
)

cryosphere = DataFrame(
    "Process" => [
        "formation of permafrost",
        "formation of sea ice",
        "formation of land ice masses"
    ],
    "Characteristic Time Scale" => [
        [1e7, 1e9],
        [1e7, 1e8],
        [1e8, 1e11],
    ],
    "Characteristic Spatial Scale" => [
        [1e0, 1e6],
        [1e0, 1e6],
        [1e2, 1e7]
    ]
)

land_surface = DataFrame(
    "Process" => [
        "changes in reflectivity",
        "isostatic equilibration of the crust by covering ice masses"
    ],
    "Characteristic Time Scale" => [
        [1e7, 1e8],
        [1e8, 1e11]
    ],
    "Characteristic Spatial Scale" => [
        [1e2, Inf],
        [1e6, Inf]
    ]
)

biosphere = DataFrame(
    "Process" => [
        "exchange of carbon with the atmosphere",
        "transformation of vegetation zones"
    ],
    "Characteristic Time Scale" => [
        [1e4, 1e8],
        [1e9, 1e10]
    ],
    "Characteristic Spatial Scale" => [
        [1e-3, Inf],
        [1e2, 1e7]
    ]
)