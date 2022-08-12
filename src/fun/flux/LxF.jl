# Rusanov numerical fluxes
@views function fluxLxF(UL,UR,FL,FR,g,c,type)
    fLxF  = zeros(Float64,3)
    for dim ∈ 1:3
            fLxF[dim] = (0.5*(FL[dim]+FR[dim])-0.5*c*(UR[dim]-UL[dim]))
    end
    return fLxF 
end
@views function LxF!(FLR,UL,UR,FL,FR,g,c,nx,ny,type)
    for j ∈ 1:ny
        for i ∈ 1:nx+1
            FLR[i,j,:] = fluxLxF(UL[i,j,:],UR[i,j,:],FL[i,j,:],FR[i,j,:],g,c,type)
        end
    end
end