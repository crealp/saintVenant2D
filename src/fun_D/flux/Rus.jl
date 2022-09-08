# Rusanov numerical fluxes
@views function Rus_D(F,UL,UR,FL,FR,g,nx,ny,type)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if i<=nx+1 && j<=ny
            hL    = max(0.0,UL[i,j,1])
            if hL>0.0
                if type == 1#"F"
                    uL = UL[i,j,2]/hL
                elseif type == 2#"G"
                    uL = UL[i,j,3]/hL
                end 
            else
                uL = 0.0
            end
            hR    = max(0.0,UR[i,j,1])
            if hR>0.0
                if type == 1#"F"
                    uR = UR[i,j,2]/hR
                elseif type == 2#"G"
                    uR = UR[i,j,3]/hR
                end 
            else
                uR = 0.0
            end
            cL    = sqrt(g*hL)
            cR    = sqrt(g*hR)
            λ     = max(abs(uL-sqrt(g*hL)),abs(uR+sqrt(g*hR)))
            for dim ∈ 1:3
                F[i,j,dim] = (0.5*(FL[i,j,dim]+FR[i,j,dim])-0.5*λ*(UR[i,j,dim]-UL[i,j,dim]))
            end
    end
    return nothing
end