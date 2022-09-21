module saintVenant
export greet
export geoflow,runoff

greet()   = print("Welcome in saintVenant module!")
geoflow() = include("./scripts/geoflow.jl")
runoff()  = include("./scripts/runoff.jl")

println(greet())
end # module
