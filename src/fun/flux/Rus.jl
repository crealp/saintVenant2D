# Rusanov numerical fluxes
@views function fluxRus(UL,UR,FL,FR,g,type)
    hL    = max(0.0,UL[1])
    if hL>0.0
        if type == "F"
            uL = UL[2]/hL
        elseif type == "G"
            uL = UL[3]/hL
        end 
    else
        uL = 0.0
    end
    hR    = max(0.0,UR[1])
    if hR>0.0
        if type == "F"
            uR = UR[2]/hR
        elseif type == "G"
            uR = UR[3]/hR
        end 
    else
        uR = 0.0
    end
    cL    = sqrt(g*hL)
    cR    = sqrt(g*hR)
    λ     = maximum(abs.(hcat(uL-sqrt(g*hL),uR+sqrt(g*hR))))
    fRus  = zeros(Float64,3)
    for dim in 1:3
        fRus[dim] = (0.5*(FL[dim]+FR[dim])-0.5*λ*(UR[dim]-UL[dim]))
    end
    return fRus 
end
@views function Rus!(FLR,UL,UR,FL,FR,g,nx,ny,type)
    for j in 1:ny
        for i in 1:nx+1
            FLR[i,j,:] = fluxRus(UL[i,j,:],UR[i,j,:],FL[i,j,:],FR[i,j,:],g,type)
        end
    end
end