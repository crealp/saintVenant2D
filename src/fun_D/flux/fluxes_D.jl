# Hydrostatic reconstruction for variable topography
@views function fluxRus_D(UFS,U,z,g,nx,ny,dim) # see well-balanced scheme, e.g., http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
    # index initialization
    i = (blockIdx().x-1)*blockDim().x+threadIdx().x
    j = (blockIdx().y-1)*blockDim().y+threadIdx().y
    if dim == 1 # flux function vector F
        if i<=nx && j<=ny
            # well-balanced scheme
            hL = max(0.0,U[i,j,1]+z[i,j]-max(z[i,j],z[i+1,j]))
            if hL > 0.0
                UFS[i,j,1,1] = hL
                UFS[i,j,2,1] = hL*U[i,j,2]/U[i,j,1]
                UFS[i,j,3,1] = hL*U[i,j,3]/U[i,j,1]
                UFS[i,j,1,3] = UFS[i,j,2,1]
                UFS[i,j,2,3] = hL*(UFS[i,j,2,1]/UFS[i,j,1,1])^2+0.5*g*UFS[i,j,1,1]^2
                UFS[i,j,3,3] = hL*(UFS[i,j,2,1]/UFS[i,j,1,1])*(UFS[i,j,3,1]/UFS[i,j,1,1])
            end
            hR  = max(0.0,U[i+1,j,1]+z[i+1,j]-max(z[i,j],z[i+1,j]))
            if hR > 0.0
                UFS[i,j,1,2] = hR 
                UFS[i,j,2,2] = hR*U[i+1,j,2]/U[i+1,j,1]
                UFS[i,j,3,2] = hR*U[i+1,j,3]/U[i+1,j,1]
                UFS[i,j,1,4] = UFS[i,j,2,2]
                UFS[i,j,2,4] = hR*(UFS[i,j,2,2]/UFS[i,j,1,2])^2+0.5*g*UFS[i,j,1,2]^2
                UFS[i,j,3,4] = hR*(UFS[i,j,2,2]/UFS[i,j,1,2])*(UFS[i,j,3,2]/UFS[i,j,1,2])
            end
            UFS[i,j,2,5] = 0.5*g*(U[i+1,j,1]^2-hR^2)
            UFS[i,j,2,6] = 0.5*g*(U[i  ,j,1]^2-hL^2)
        end
        # Rusanov flux definition
        if i<=nx+1 && j<=ny
            hL = max(0.0,UFS[i,j,1,1])
            if hL>0.0
                uL = UFS[i,j,2,1]/hL
            else
                uL = 0.0
            end
            hR = max(0.0,UFS[i,j,1,2])
            if hR>0.0
                uR = UFS[i,j,2,2]/hR
            else
                uR = 0.0
            end
            cL = sqrt(g*hL)
            cR = sqrt(g*hR)
            λ  = max(abs(uL-sqrt(g*hL)),abs(uR+sqrt(g*hR)))
            for dim ∈ 1:3
                UFS[i,j,dim,7] = (0.5*(UFS[i,j,dim,3]+UFS[i,j,dim,4])-0.5*λ*(UFS[i,j,dim,2]-UFS[i,j,dim,1]))
            end
        end
    elseif dim == 2 # flux function vector G
        if i<=nx && j<=ny
            hL = max(0.0,U[i,j,1]+z[i,j]-max(z[i,j],z[i,j+1]))
            if hL > 0.0
                UFS[i,j,1,1] = hL
                UFS[i,j,2,1] = hL*U[i,j,2]/U[i,j,1]
                UFS[i,j,3,1] = hL*U[i,j,3]/U[i,j,1]
                UFS[i,j,1,3] = UFS[i,j,3,1]
                UFS[i,j,2,3] = hL*(UFS[i,j,2,1]/UFS[i,j,1,1])*(UFS[i,j,3]/UFS[i,j,1,1])
                UFS[i,j,3,3] = hL*(UFS[i,j,3,1]/UFS[i,j,1,1])^2+0.5*g*UFS[i,j,1,1]^2
            end
            hR  = max(0.0,U[i,j+1,1]+z[i,j+1]-max(z[i,j],z[i,j+1]))
            if hR > 0.0
                UFS[i,j,1,2] = hR 
                UFS[i,j,2,2] = hR*U[i,j+1,2]/U[i,j+1,1]
                UFS[i,j,3,2] = hR*U[i,j+1,3]/U[i,j+1,1]
                UFS[i,j,1,4] = UFS[i,j,3,2]
                UFS[i,j,2,4] = hR*(UFS[i,j,2,2]/UFS[i,j,1,2])*(UFS[i,j,3,2]/UFS[i,j,1,2])
                UFS[i,j,3,4] = hR*(UFS[i,j,3,2]/UFS[i,j,1,2])^2+0.5*g*UFS[i,j,1,2]^2
            end
            UFS[i,j,3,5] = 0.5*g*(U[i,j+1,1]^2-hR^2)
            UFS[i,j,3,6] = 0.5*g*(U[i  ,j,1]^2-hL^2)
        end
        # Rusanov flux definition
        if i<=nx && j<=ny+1
            hL = max(0.0,UFS[i,j,1,1])
            if hL>0.0
                uL = UFS[i,j,3,1]/hL
            else
                uL = 0.0
            end
            hR = max(0.0,UFS[i,j,1,2])
            if hR>0.0
                uR = UFS[i,j,3,2]/hR
            else
                uR = 0.0
            end
            cL = sqrt(g*hL)
            cR = sqrt(g*hR)
            λ  = max(abs(uL-sqrt(g*hL)),abs(uR+sqrt(g*hR)))
            for dim ∈ 1:3
                UFS[i,j,dim,7] = (0.5*(UFS[i,j,dim,3]+UFS[i,j,dim,4])-0.5*λ*(UFS[i,j,dim,2]-UFS[i,j,dim,1]))
            end
        end
    end    
    return nothing
end