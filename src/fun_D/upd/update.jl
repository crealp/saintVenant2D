@views function updateU_D(U,F,SL,SR,c,nx,ny,dim)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if i<=nx && j<=ny && dim==1
        U[i,j,1] -= c*((F[i+1,j  ,1].+SR[i+1,j  ,1]).-(F[i,j,1].+SL[i,j,1]))
        U[i,j,2] -= c*((F[i+1,j  ,2].+SR[i+1,j  ,2]).-(F[i,j,2].+SL[i,j,2]))
        U[i,j,3] -= c*((F[i+1,j  ,3].+SR[i+1,j  ,3]).-(F[i,j,3].+SL[i,j,3]))
    elseif i<=nx && j<=ny && dim==2
        U[i,j,1] -= c*((F[i  ,j+1,1].+SR[i  ,j+1,1]).-(F[i,j,1].+SL[i,j,1]))
        U[i,j,2] -= c*((F[i  ,j+1,2].+SR[i  ,j+1,2]).-(F[i,j,2].+SL[i,j,2]))
        U[i,j,3] -= c*((F[i  ,j+1,3].+SR[i  ,j+1,3]).-(F[i,j,3].+SL[i,j,3]))
    end
    return nothing
end
@views function updateAdvU_D(U,S,Δt,nx,ny)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if i<=nx && j<=ny
        U[i,j,1]+=Δt*S[i,j,1]
        U[i,j,2]+=Δt*S[i,j,2]
        U[i,j,3]+=Δt*S[i,j,3]
    end
    return nothing
end