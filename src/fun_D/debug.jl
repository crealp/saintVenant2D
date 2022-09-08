# copy-paste in julia REPL

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
    tC     = 1.0./25.0
    type   = "HLL"
    U,F,G  = getUF(h,Qx,Qy,g,nx,ny)
    dim    = 3
    Δt     = 0.1


    @code_warntype fluxes(U,F,z,g,type,nx,ny,dim)
    @code_warntype wellBal(U,z,g,nx,ny,3,"F")
    @code_warntype advSolve(h,Qx,Qy,z,U,F,G,g,Δx,Δy,Δt,nx,ny,type)