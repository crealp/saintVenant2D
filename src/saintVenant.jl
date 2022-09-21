module saintVenant
export geoflow,runoff,coast,basin


# include dependencies & function call(s)
include(joinpath("./fun", "superInclude.jl"))

# include geoflow routine in saintVenant module
@doc raw"""
    geoflow(lx::Float64,ly::Float64,nx::Int64): solves a non-linear hyperbolic 2D Saint-Venant problem considering a Coulomb-type rheology within a finite volume framework on a Cartesian grid
    # args:
    - lx : dimension along the x-direciton.
    - ly : dimension along the y-direciton.
    - nx : number of grid nodes along the x-direction.
"""
geoflow()
include(joinpath("../scripts", "geoflow.jl"))


include(joinpath("../scripts", "runoff.jl"))
include(joinpath("../scripts", "coast.jl"))
include(joinpath("../scripts", "basin.jl"))

end # module