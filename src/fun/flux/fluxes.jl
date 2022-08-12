include("wellBal.jl")
include("HLL.jl")
include("HLLC.jl")
include("Rus.jl")
include("LxF.jl")
# Numerical fluxes function
@views function fluxes!(F,U,z,g,Δx,Δt,type,nx,ny,dim)
    # out 
    f  = zeros(Float64,nx+1,ny,3)
    sl = zeros(Float64,nx+1,ny,3)
    sr = zeros(Float64,nx+1,ny,3)
    # in
    UL = zeros(Float64,nx+1,ny,3)
    UR = zeros(Float64,nx+1,ny,3)
    FL = zeros(Float64,nx+1,ny,3)
    FR = zeros(Float64,nx+1,ny,3)    
    SL = zeros(Float64,nx+1,ny,1)
    SR = zeros(Float64,nx+1,ny,1)

    if dim=="x"
        wellBal!(UL,UR,FL,FR,SL,SR,U,z,g,nx,ny,3,"F")
        sl[:,:,2].=SL[:,:]
        sr[:,:,2].=SR[:,:]
    elseif dim=="y"
        wellBal!(UL,UR,FL,FR,SL,SR,U,z,g,nx,ny,3,"G")
        sl[:,:,3].=SL[:,:]
        sr[:,:,3].=SR[:,:]
    end
    if type=="Rus"
        if dim=="x"
            Rus!(f,UL,UR,FL,FR,g,nx,ny,"F")
        elseif dim=="y"
            Rus!(f,UL,UR,FL,FR,g,nx,ny,"G")
        end
    elseif type=="HLL"
        if dim=="x"
            HLL!(f,UL,UR,FL,FR,g,nx,ny,"F")
        elseif dim=="y"
            HLL!(f,UL,UR,FL,FR,g,nx,ny,"G")
        end
    elseif type=="HLLC"
        if dim=="x"
            HLLC!(f,UL,UR,FL,FR,g,nx,ny,"F")
        elseif dim=="y"
            HLLC!(f,UL,UR,FL,FR,g,nx,ny,"G")
        end    
    else 
        println("no numerical fluxes defined !")
        println("available flux functions are: a) Rus  - Rusanov fluxes")
        println("                            : b) HLL  - HLL approximate Riemann solver")
        println("                            : c) HLLC - HLLC  approximate Riemann solver")
        exit(-1)
    end
    F.=(f[2:nx+1,:,:].+sr[2:nx+1,:,:]).-(f[1:nx+0,:,:].+sl[1:nx+0,:,:])
end