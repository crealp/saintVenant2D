module saintVenant
export greet
export geoflow,runoff,coast,basin

greet()   = print("Welcome in saintVenant module!")
geoflow() = include("./scripts/geoflow.jl")
runoff()  = include("./scripts/runoff.jl")
coast()   = include("./scripts/coast.jl")
basin()   = include("./scripts/basin.jl")

end # module

#geoflow() = main()
#include(joinpath("../scripts", "geoflow.jl"))