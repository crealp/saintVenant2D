# initialisation & global definition(s)
using CUDA
# include dependencies
	include(joinpath("./bc"   , "getBCs_D.jl"))
	include(joinpath("./flux" , "fluxes_D.jl"))
	include(joinpath("./solve", "advSolve_D.jl"))
	include(joinpath("./solve", "get_D.jl"))
	include(joinpath("./solve", "souSolve_D.jl"))
	include(joinpath("./solve", "svSolver_D.jl"))
	include(joinpath("./upd"  , "update_D.jl"))
