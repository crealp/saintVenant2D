@views function getUF(h,Qx,Qy,g,nx,ny)
    U = zeros(Float64,nx,ny,3)
    F = zeros(Float64,nx,ny,3)
    G = zeros(Float64,ny,nx,3)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j] > 0.0
                U[i,j,1] = h[i,j]
                U[i,j,2] = Qx[i,j]
                U[i,j,3] = Qy[i,j]

                F[i,j,1] = Qx[i,j]
                F[i,j,2] = h[i,j]*(Qx[i,j]/h[i,j])^2+0.5*g*h[i,j]^2
                F[i,j,3] = h[i,j]*(Qx[i,j]/h[i,j])*(Qy[i,j]/h[i,j])

                G[j,i,1] = Qy[i,j]
                G[j,i,2] = h[i,j]*(Qx[i,j]/h[i,j])*(Qy[i,j]/h[i,j])
                G[j,i,3] = h[i,j]*(Qy[i,j]/h[i,j])^2+0.5*g*h[i,j]^2
            end
        end
    end
    return U,F,G
end
@views function getU(h,Qx,Qy,nx,ny)
    U = zeros(Float64,nx,ny,3)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j] > 0.0
                U[i,j,1] = h[i,j]
                U[i,j,2] = Qx[i,j]
                U[i,j,3] = Qy[i,j]
            end
        end
    end
    return U
end
@views function getU!(U,h,Qx,Qy,nx,ny)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j] > 0.0
                U[i,j,1] = h[i,j]
                U[i,j,2] = Qx[i,j]
                U[i,j,3] = Qy[i,j]
            else
                U[i,j,1] = 0.0
                U[i,j,2] = 0.0
                U[i,j,3] = 0.0
            end
        end
    end
    return nothing
end
@views function getΔt(h,Qx,Qy,g,Δx,Δy,CFL,nx,ny)
    # find minimal Δt consistent with CFL
    cx = 0.0
    cy = 0.0
    Δ  = min(Δx,Δy)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j] > 0.0
                u  = Qx[i,j]/(h[i,j])
                v  = Qy[i,j]/(h[i,j])
                cx = max(cx,abs(u)+sqrt(g*h[i,j]))
                cy = max(cy,abs(v)+sqrt(g*h[i,j]))
            end
        end
    end
    return CFL*Δ/(max(cx,cy))
end