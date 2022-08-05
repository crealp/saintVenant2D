# include dependencies & function call(s)
include("../src/fun/usingPackages.jl")
include("../src/fun/misc.jl")
include("../src/fun/geometry.jl")
include("../src/fun/solve/svSolver.jl")

@views function main()
    Dsim   = param("HLLC",
                    false,
                    "newtonian",
                    true
                )
    # physical constant
    g     = 9.81
    # number of points
    path = "./docs/example/hillshade/data/dtm_1m.txt"
    z = CSV.read(path,DataFrame,header=false; delim="\t")
    show(z)
    z = Array(z)'
    nx,ny = size(z)
    z = z[200:round(Int64,nx/2),round(Int64,ny/2):end]


    nx,ny = size(z)
    Δx    = 1.0
    Δy    = Δx
    xc    = 1:Δx:nx*Δx
    yc    = 1:Δx:ny*Δy
    h     = 1.0e-3.*ones(Float64,nx,ny)

    Qx    = zeros(Float64,nx,ny)
    Qy    = zeros(Float64,nx,ny)
    # action
    CFL   = 0.5
    T     = 30.0*60.0
    tC    = 1.0
    svSolver(xc,yc,h,Qx,Qy,z,g,CFL,T,tC,Δx,Δy,nx,ny,Dsim)
end
main()
# https://techytok.com/lesson-parallel-computing
# https://nbviewer.org/github/daniel-koehn/Differential-equations-earth-system/blob/master/10_Shallow_Water_Equation_2D/01_2D_Shallow_Water_Equations.ipynb