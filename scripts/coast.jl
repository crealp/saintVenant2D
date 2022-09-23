@views function coast()
    Dsim   = param("HLLC",
                    false,
                    "newtonian",
                    false
                )
    # physical constant
    lx     = 20.0
    ly     = 10.0
    g      = 9.81
    # number of points
    nx     = 200
    ny     = Int64((ly/lx)*nx)
    Qx     = zeros(Float64,nx,ny)
    Qy     = zeros(Float64,nx,ny)
    h,z,xc,yc,Δx,Δy = hyperbolic_floor(lx,ly,nx,ny)
    # action
    CFL = 0.5
    T  = 12.0
    tC = 1.0./25.0
    svSolver(xc,yc,h,Qx,Qy,z,g,CFL,T,tC,Δx,Δy,nx,ny,Dsim)
end