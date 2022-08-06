# Hydrostatic reconstruction for variable topography
@views function wellBal!(UL,UR,FL,FR,SL,SR,U,z,g,nx,ny,K,dim) # see well-balanced scheme, e.g., http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
    if dim=="F"
        for j in 1:ny
            for i in 1:nx+1
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
            SR[i,j] = 0.5*g*(U[i  ,j,1]^2-hL^2)
            SL[i,j] = 0.5*g*(U[i+1,j,1]^2-hR^2)
            end
        end
    elseif dim=="G"
        for j in 1:ny
            for i in 1:nx+1
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
            SR[i,j] = 0.5*g*(U[i  ,j,1]^2-hL^2)
            SL[i,j] = 0.5*g*(U[i+1,j,1]^2-hR^2)
            end
        end
    end    
end
@views function wellBal2!(UL,UR,FL,FR,SL,SR,U,z,g,nx,ny,K,dim) # see well-balanced scheme, e.g., http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
    ϕb = 10.0*π/180
    ϕi = 35.0*π/180
    Kact = 2*(1.0-sqrt(1.0-cos(ϕi)^2*(1.0+tan(ϕb)^2)))/cos(ϕi)^2-1
    Kpas = 2*(1.0+sqrt(1.0-cos(ϕi)^2*(1.0+tan(ϕb)^2)))/cos(ϕi)^2-1
    for j in 1:ny
        for i in 1:nx+1
        hL = max(0.0,U[i,j,1]+z[i,j]-max(z[i,j],z[i+1,j]))
        if hL > 0.0
            UL[i,j,1] = hL
            UL[i,j,2] = hL*U[i,j,2]/U[i,j,1]
            UL[i,j,3] = hL*U[i,j,3]/U[i,j,1]
        end
        hR  = max(0.0,U[i+1,j,1]+z[i+1,j]-max(z[i,j],z[i+1,j]))
        if hR > 0.0
            UR[i,j,1] = hR 
            UR[i,j,2] = hR*U[i+1,j,2]/U[i+1,j,1]
            UR[i,j,3] = hR*U[i+1,j,3]/U[i+1,j,1]
        end

        if (UR[i,j,2]-UL[i,j,2])>0.0
            kx = Kact
        elseif (UR[i,j,2]-UL[i,j,2])<0.0
            kx = Kpas
        else
            kx = 1.0
        end
        if (UR[i,j,3]-UL[i,j,3])>0.0
            ky = Kact
        elseif (UR[i,j,3]-UL[i,j,3])<0.0
            ky = Kpas
        else
            ky = 1.0
        end

        if hL > 0.0
            if dim == "F"
                FL[i,j,1] = UL[i,j,2]
                FL[i,j,2] = hL*(UL[i,j,2]/UL[i,j,1])^2+kx*0.5*g*UL[i,j,1]^2
                FL[i,j,3] = hL*(UL[i,j,2]/UL[i,j,1])*(UL[i,j,3]/UL[i,j,1])
            elseif dim == "G"
                FL[i,j,1] = UL[i,j,3]
                FL[i,j,2] = hL*(UL[i,j,2]/UL[i,j,1])*(UL[i,j,3]/UL[i,j,1])
                FL[i,j,3] = hL*(UL[i,j,3]/UL[i,j,1])^2+ky*0.5*g*UL[i,j,1]^2
            end
        end
        if hR > 0.0
            if dim == "F"
                FR[i,j,1] = UR[i,j,2]
                FR[i,j,2] = hR*(UR[i,j,2]/UR[i,j,1])^2+kx*0.5*g*UR[i,j,1]^2
                FR[i,j,3] = hR*(UR[i,j,2]/UR[i,j,1])*(UR[i,j,3]/UR[i,j,1])
            elseif dim == "G"
                FR[i,j,1] = UR[i,j,3]
                FR[i,j,2] = hR*(UR[i,j,2]/UR[i,j,1])*(UR[i,j,3]/UR[i,j,1])
                FR[i,j,3] = hR*(UR[i,j,3]/UR[i,j,1])^2+ky*0.5*g*UR[i,j,1]^2
            end
        end

        SR[i,j] = 0.5*g*(U[i  ,j,1]^2-hL^2)
        SL[i,j] = 0.5*g*(U[i+1,j,1]^2-hR^2)
        end
    end
end