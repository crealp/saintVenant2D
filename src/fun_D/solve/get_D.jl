@views function getQxQyh_D(h,Qx,Qy,U,g,nx,ny)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    # calculation
    if i<=nx && j<=ny
        h[i,j] = U[i,j,1] 
        Qx[i,j]= U[i,j,2]
        Qy[i,j]= U[i,j,3]
    end
    return nothing
end
@views function getUF_D(U,F,G,h,Qx,Qy,g,nx,ny)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    # calculation
    if i<=nx && j<=ny
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
    return nothing
end
@views function getU_D(U,h,Qx,Qy,nx,ny)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    # calculation
    if i<=nx && j<=ny
        if h[i,j] > 0.0
            U[i,j,1] = h[i,j]
            U[i,j,2] = Qx[i,j]
            U[i,j,3] = Qy[i,j]
        end
    end
    return nothing
end
@views function setUFS_D(UFS,nx,ny)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    # calculation
    if i<=nx && j<=ny
        k=1
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
        k=2
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
        k=3
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
        k=4
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
        k=5
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
        k=6
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
        k=7
            UFS[i,j,1,k] = 0.0
            UFS[i,j,2,k] = 0.0
            UFS[i,j,3,k] = 0.0
    end
    return nothing
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