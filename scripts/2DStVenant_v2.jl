# include dependencies & function call(s)
include("../src/fun/usingPackages.jl")
include("../src/fun/misc.jl")
include("../src/fun/geometry.jl")
include("../src/fun/solve/svSolver.jl")

@views function main()
    Dsim   = param("HLLC",false,"coulomb",false)
    #Dsim   = param("HLLC",false,"newtonian",false)
    # physical constant
    lx     = 20.0
    ly     = 20.0
    g      = 9.81
    # number of points
    nx     = 200
    ny     = Int64((ly/lx)*nx)
    Qx     = zeros(Float64,nx,ny)
    Qy     = zeros(Float64,nx,ny)
    #h,z,xc,yc,Δx,Δy = bowl_floor_noΔz(lx,ly,nx,ny)
    h,z,xc,yc,Δx,Δy = bowl_floor(lx,ly,nx,ny)
    # action
    CFL    = 0.5
    T      = 5.0
    tC     = 1.0./25
    svSolver(xc,yc,h,Qx,Qy,z,g,CFL,T,tC,Δx,Δy,nx,ny,Dsim)
end
main()
# https://techytok.com/lesson-parallel-computing
# https://nbviewer.org/github/daniel-koehn/Differential-equations-earth-system/blob/master/10_Shallow_Water_Equation_2D/01_2D_Shallow_Water_Equations.ipynb