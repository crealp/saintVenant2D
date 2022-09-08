# Hydrostatic reconstruction for variable topography
@views function wellBal_D(UL,UR,FL,FR,SL,SR,U,z,g,nx,ny,dim) # see well-balanced scheme, e.g., http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if dim == 1#"F"
        if i<=nx && j<=ny
            hL = max(0.0,U[i,j,1]+z[i,j]-max(z[i,j],z[i+1,j]))
            if hL > 0.0
                UL[i,j,1] = hL
                UL[i,j,2] = hL*U[i,j,2]/U[i,j,1]
                UL[i,j,3] = hL*U[i,j,3]/U[i,j,1]
                FL[i,j,1] = UL[i,j,2]
                FL[i,j,2] = hL*(UL[i,j,2]/UL[i,j,1])^2+0.5*g*UL[i,j,1]^2
                FL[i,j,3] = hL*(UL[i,j,2]/UL[i,j,1])*(UL[i,j,3]/UL[i,j,1])
            end
            hR  = max(0.0,U[i+1,j,1]+z[i+1,j]-max(z[i,j],z[i+1,j]))
            if hR > 0.0
                UR[i,j,1] = hR 
                UR[i,j,2] = hR*U[i+1,j,2]/U[i+1,j,1]
                UR[i,j,3] = hR*U[i+1,j,3]/U[i+1,j,1]
                FR[i,j,1] = UR[i,j,2]
                FR[i,j,2] = hR*(UR[i,j,2]/UR[i,j,1])^2+0.5*g*UR[i,j,1]^2
                FR[i,j,3] = hR*(UR[i,j,2]/UR[i,j,1])*(UR[i,j,3]/UR[i,j,1])
            end
            SR[i,j,2] = 0.5*g*(U[i  ,j,1]^2-hL^2)
            SL[i,j,2] = 0.5*g*(U[i+1,j,1]^2-hR^2)
        end
    elseif dim == 2#"G"
        if i<=nx && j<=ny
            hL = max(0.0,U[i,j,1]+z[i,j]-max(z[i,j],z[i+1,j]))
            if hL > 0.0
                UL[i,j,1] = hL
                UL[i,j,2] = hL*U[i,j,2]/U[i,j,1]
                UL[i,j,3] = hL*U[i,j,3]/U[i,j,1]
                FL[i,j,1] = UL[i,j,3]
                FL[i,j,2] = hL*(UL[i,j,2]/UL[i,j,1])*(UL[i,j,3]/UL[i,j,1])
                FL[i,j,3] = hL*(UL[i,j,3]/UL[i,j,1])^2+0.5*g*UL[i,j,1]^2
            end
            hR  = max(0.0,U[i+1,j,1]+z[i+1,j]-max(z[i,j],z[i+1,j]))
            if hR > 0.0
                UR[i,j,1] = hR 
                UR[i,j,2] = hR*U[i+1,j,2]/U[i+1,j,1]
                UR[i,j,3] = hR*U[i+1,j,3]/U[i+1,j,1]
                FR[i,j,1] = UR[i,j,3]
                FR[i,j,2] = hR*(UR[i,j,2]/UR[i,j,1])*(UR[i,j,3]/UR[i,j,1])
                FR[i,j,3] = hR*(UR[i,j,3]/UR[i,j,1])^2+0.5*g*UR[i,j,1]^2
            end
            SR[i,j,3] = 0.5*g*(U[i  ,j,1]^2-hL^2)
            SL[i,j,3] = 0.5*g*(U[i+1,j,1]^2-hR^2)
        end
    end    
    return nothing
end