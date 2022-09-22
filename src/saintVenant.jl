module saintVenant
export geoflow,runoff,coast,basin
# include dependencies & function call(s)
include(joinpath("./fun", "superInclude.jl"))
global path_plot = "viz/out/"
    if isdir(path_plot)==false
        mkdir(path_plot)    
    end
global path_save = "viz/dat/"
    if isdir(path_save)==false
        mkdir(path_save)    
    end
# include geoflow routine in saintVenant module
@doc raw"""
    geoflow(lx::Float64,ly::Float64,nx::Int64,rheoType::String,solvType::String): solves a non-linear hyperbolic 2D Saint-Venant problem considering a Coulomb-type rheology within a finite volume framework on a Cartesian grid
    # args:
    - lx       : dimension along the x-direciton.
    - ly       : dimension along the y-direciton.
    - nx       : number of grid nodes along the x-direction.
    - rheoType : select the rheology, i.e., "coulomb", "newtonian" or "plastic"
    - solveType: select the numerical flux, i.e., "Rusanov", "HLL" or "HLLC"
"""
geoflow()
include(joinpath("../scripts", "geoflow.jl"))
# include runoff routine in saintVenant module
include(joinpath("../scripts", "runoff.jl"))
# include coast routine in saintVenant module
include(joinpath("../scripts", "coast.jl"))
# include basin routine in saintVenant module
include(joinpath("../scripts", "basin.jl"))
end

#=
----------------------------------------------------------------------
                **FOR DEVELOPMENT**

When changing stuffs within the module, in REPL, enter the following:
    julia> include("./src/saintVenant.jl")
    WARNING: replacing module saintVenant.
    Main.saintVenant
    julia> saintVenant.geoflow(20.0,10.0,200)

----------------------------------------------------------------------
=#