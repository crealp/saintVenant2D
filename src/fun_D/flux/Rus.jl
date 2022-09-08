# Rusanov numerical fluxes
@views function Rus_D(UFS,g,nx,ny,type)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if i<=nx+1 && j<=ny
            hL    = max(0.0,UFS[i,j,1,1])
            if hL>0.0
                if type == 1#"F"
                    uL = UFS[i,j,2,1]/hL
                elseif type == 2#"G"
                    uL = UFS[i,j,3,1]/hL
                end 
            else
                uL = 0.0
            end
            hR    = max(0.0,UFS[i,j,1,2])
            if hR>0.0
                if type == 1#"F"
                    uR = UFS[i,j,2,2]/hR
                elseif type == 2#"G"
                    uR = UFS[i,j,3,2]/hR
                end 
            else
                uR = 0.0
            end
            cL    = sqrt(g*hL)
            cR    = sqrt(g*hR)
            λ     = max(abs(uL-sqrt(g*hL)),abs(uR+sqrt(g*hR)))
            for dim ∈ 1:3
                UFS[i,j,dim,7] = (0.5*(UFS[i,j,dim,3]+UFS[i,j,dim,4])-0.5*λ*(UFS[i,j,dim,2]-UFS[i,j,dim,1]))
            end
    end
    return nothing
end