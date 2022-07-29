# approximate HLLC Riemann solver, see Creed etal, 2016
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
        uL = 0.0
        cR = sqrt(g*hR)
        sL = uR-2.0*cR
        sR = uR+1.0*cR
    elseif hR==0.0 && hL>0.0 # dry-bed case on the right  
        if type == "F"
            uL = UL[2]/hL
        elseif type == "G"
            uL = UL[3]/hL
        end 
        uR = 0.0
        cL = sqrt(g*hL)
        sL = uL-1.0*cL
        sR = uL+2.0*cL 
    else
        uL = 0.0
        uR = 0.0
        sL = 0.0
        sR = 0.0
    end   
    sS = (sL*hR*(uR-sR)-sR*hL*(uL-sL))/(hR*(uR-sR)-hL*(uL-sL))

    vL = 0.0
    vR = 0.0
    if hL>0.0
        if type == "F"
            vL = UL[3]/hL
        elseif type == "G"
            vL = UL[2]/hL
        end 
    end
    if hR>0.0
        if type == "F"
            vR = UR[3]/hR
        elseif type == "G"
            vR = UR[2]/hR
        end 
    end

    return sL,sR,sS,vL,vR 
end
@views function fluxHLLC(fHLL,fstarL,fstarR,fstar,UL,UR,FL,FR,g,type,nD)
    sL,sR,sS,vL,vR = sLsR(UL,UR,g,type)
    fstar         .= (sR.*FL.-sL.*FR.+sL.*sR.*(UR.-UL))./(sR.-sL)
    if type == "F"
        fstarL[1]=fstar[1]
        fstarL[2]=fstar[2]
        fstarL[3]=fstar[1]*vL
        
        fstarR[1]=fstar[1]
        fstarR[2]=fstar[2]
        fstarR[3]=fstar[1]*vR
    elseif type == "G"
        fstarL[1]=fstar[1]
        fstarL[2]=fstar[1]*vL
        fstarL[3]=fstar[3]

        fstarR[1]=fstar[1]
        fstarR[2]=fstar[1]*vR
        fstarR[3]=fstar[3]
    end
    for dim ∈ 1:nD
        if     sL>=0.0
            fHLL[dim] = FL[dim]
        elseif sL<0.0<=sS
            fHLL[dim] = fstarL[dim]
        elseif sS<0.0<=sR
            fHLL[dim] = fstarR[dim]
        elseif sR<=0.0
            fHLL[dim] = FR[dim]
        end
    end
    return fHLL
end
@views function HLLC!(F,UL,UR,FL,FR,g,nx,ny,type)
    fHLL   = zeros(Float64,3)
    fstarL = zeros(Float64,3)
    fstarR = zeros(Float64,3)
    fstar  = zeros(Float64,3)
    for j ∈ 1:ny
        for i ∈ 1:nx+1
            F[i,j,:] .= fluxHLLC(fHLL,fstarL,fstarR,fstar,UL[i,j,:],UR[i,j,:],FL[i,j,:],FR[i,j,:],g,type,3)
        end
    end
end