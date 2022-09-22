@views function geoflow(lx::Float64,ly::Float64,nx::Int64,rheoType::String,solvType::String)
    # Dsim definition
    Dsim   = param(solvType,false,rheoType,false)
    #Dsim   = param("HLLC",false,"newtonian",false)
    # physical constant
    g      = 9.81
    # number of points
    ny     = Int64((ly/lx)*nx)
    Qx     = zeros(Float64,nx,ny)
    Qy     = zeros(Float64,nx,ny)
    #h,z,xc,yc,Δx,Δy = bowl_floor_noΔz(lx,ly,nx,ny)
    #h,z,xc,yc,Δx,Δy = bowl_floor(lx,ly,nx,ny)
    h,z,xc,yc,Δx,Δy = incline(lx,ly,nx,ny)
    # action
    CFL    = 0.5
    T      = 5.0
    tC     = 1.0/25.0
    svSolver(xc,yc,h,Qx,Qy,z,g,CFL,T,tC,Δx,Δy,nx,ny,Dsim)
end