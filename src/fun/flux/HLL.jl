# approximate HLL Riemann solver, see Castro-Orgaz etal, 2019, Shallow Water Hydraulics, pp. 363-365
@views function sLsR(UL,UR,g,type)
    hL = max(0.0,UL[1])
    hR = max(0.0,UR[1])
    if hL>0.0 && hR>0.0
        if type == "F"
            uL = UL[2]/hL
            uR = UR[2]/hR
        elseif type == "G"
            uL = UL[3]/hL
            uR = UR[3]/hR
        end 
        cL = sqrt(g*hL)
        cR = sqrt(g*hR)
        hS = ((0.5*(cL+cR)+0.25*(uL-uR))^2)/g
        uS = 0.5*(uL+uR)+cL-cR
        cS = sqrt(g*hS)
        sL = min(uL-cL,uS-cS)
        sR = max(uR+cR,uS+cS)
    elseif hL==0.0 && hR>0.0 # dry-bed case on the left  
        if type == "F"
            uR = UR[2]/hR
        elseif type == "G"
            uR = UR[3]/hR
        end 
        cR = sqrt(g*hR)
        sL = uR-2.0*cR
        sR = uR+1.0*cR
    elseif hR==0.0 && hL>0.0 # dry-bed case on the right  
        if type == "F"
            uL = UL[2]/hL
        elseif type == "G"
            uL = UL[3]/hL
        end 
        cL = sqrt(g*hL)
        sL = uL-1.0*cL
        sR = uL+2.0*cL 
    else
        sL = 0.0
        sR = 0.0
    end   
    return sL,sR 
end
@views function fluxHLL(HLL,UL,UR,FL,FR,g,type)
    sL,sR = sLsR(UL,UR,g,type)
    for dim ∈ 1:3
        if sL >= 0.0
            HLL[dim] = FL[dim]
        elseif sL<0.0<sR
            HLL[dim] = ( sR*FL[dim]-sL*FR[dim]+sR*sL*(UR[dim]-UL[dim]) )/(sR-sL)
        elseif sR <= 0.0
            HLL[dim] = FR[dim]
        end
    end
    return HLL
end
@views function HLL!(FLR,UL,UR,FL,FR,g,nx,ny,type)
    HLL = zeros(Float64,3)
    for j ∈ 1:ny
        for i ∈ 1:nx+1
            FLR[i,j,:].=fluxHLL(HLL,UL[i,j,:],UR[i,j,:],FL[i,j,:],FR[i,j,:],g,type)
        end
    end
    return nothing
end














































# RXIV
@views function OLD_sLsR(UL,UR,g,type)
    hL = max(0.0,UL[1])
    hR = max(0.0,UR[1])
    if hL>0.0 && hR>0.0
        if type == "F"
            uL = UL[2]/hL
            uR = UR[2]/hR
        elseif type == "G"
            uL = UL[3]/hL
            uR = UR[3]/hR
        end 
        cL = sqrt(g*hL)
        cR = sqrt(g*hR)
        hS = ((0.5*(cL+cR)+0.25*(uL-uR))^2)/g
        println(round.(vcat(uL,uR,cL,cR,hS)))
        if hS > hL
            λL = sqrt(0.5*((hS*(hS+hL))/(hL^2)))   
        elseif hS <= hL
            λL = 1.0
        end
        if hS > hR
            λR = sqrt(0.5*((hS*(hS+hR))/(hR^2)))    
        elseif hS <= hR  
            λR = 1.0   
        end
        sL = uL - cL*λL
        sR = uR + cR*λR
    elseif hL==0.0 && hR>0.0 # dry-bed case on the left  
        if type == "F"
            uR = UR[2]/hR
        elseif type == "G"
            uR = UR[3]/hR
        end 
        cR = sqrt(g*hR)
        sL = uR-2.0*cR
        sR = uR+1.0*cR
    elseif hR==0.0 && hL>0.0 # dry-bed case on the right  
        if type == "F"
            uL = UL[2]/hL
        elseif type == "G"
            uL = UL[3]/hL
        end 
        cL = sqrt(g*hL)
        sL = uL-1.0*cL
        sR = uL+2.0*cL 
    else
        sL = 0.0
        sR = 0.0
    end   
    return sL,sR 
end