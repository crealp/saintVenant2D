@views function updateU_D(U,UFS,c,nx,ny,dim)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if i<=nx && j<=ny && dim==1
        U[i,j,1] -= c*((UFS[i+1,j  ,1,7].+UFS[i+1,j  ,1,6]).-(UFS[i,j,1,7].+UFS[i,j,1,5]))
        U[i,j,2] -= c*((UFS[i+1,j  ,2,7].+UFS[i+1,j  ,2,6]).-(UFS[i,j,2,7].+UFS[i,j,2,5]))
        U[i,j,3] -= c*((UFS[i+1,j  ,3,7].+UFS[i+1,j  ,3,6]).-(UFS[i,j,3,7].+UFS[i,j,3,5]))
    elseif i<=nx && j<=ny && dim==2
        U[i,j,1] -= c*((UFS[i  ,j+1,1,7].+UFS[i  ,j+1,1,6]).-(UFS[i,j,1,7].+UFS[i,j,1,5]))
        U[i,j,2] -= c*((UFS[i  ,j+1,2,7].+UFS[i  ,j+1,2,6]).-(UFS[i,j,2,7].+UFS[i,j,2,5]))
        U[i,j,3] -= c*((UFS[i  ,j+1,3,7].+UFS[i  ,j+1,3,6]).-(UFS[i,j,3,7].+UFS[i,j,3,5]))
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