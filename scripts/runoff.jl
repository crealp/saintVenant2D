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
    path = "../dat/dtm_5m/dsm_sion.asc"
    d = Array(CSV.read(path,DataFrame,header=false; delim="\t", limit=6))
    x = Float64(d[3])
    y = Float64(d[4])
    Δ = Float64(d[5])
    z = Array(CSV.read(path,DataFrame,header=false; delim=" ", skipto=8)[:,2:end])
    z = (z[1:end,180:end])'
    nx,ny = size(z)
    xc    = 0.0:Δ:nx*Δ
    yc    = 0.0:Δ:ny*Δ

    xm  = [1750.0,2100.0]
    ym  = [750.0,1500.0]
    xm  = [0.0,1000.0]
    ym  = [0.0,1000.0]
    xf  = vcat(xm,reverse(xm))
    yf  = vcat(ym,reverse(ym))

    xId = findall(x->x>xm[1] && x<xm[2],xc)
    yId = findall(x->x>ym[1] && x<ym[2],yc)

    z0    = copy(z[xId,yId])
    xc0   = copy(xc[xId])
    yc0   = copy(yc[yId])
    Δx    = Δ
    Δy    = Δx
    nx,ny = size(z0)
    xc    = 0.0:Δx:nx*Δx
    yc    = 0.0:Δy:ny*Δy

    h     = 1.0e-6.*ones(Float64,nx,ny)
    Qx    = zeros(Float64,nx,ny)
    Qy    = zeros(Float64,nx,ny)
    # action
    CFL   = 0.5
    T     = 60.0*60.0
    tC    = 600.0
    svSolverPerf(xc0,yc0,h,Qx,Qy,z0,g,CFL,T,tC,Δx,Δy,nx,ny,Dsim)
end
main()
# https://techytok.com/lesson-parallel-computing
# https://nbviewer.org/github/daniel-koehn/Differential-equations-earth-system/blob/master/10_Shallow_Water_Equation_2D/01_2D_Shallow_Water_Equations.ipynb