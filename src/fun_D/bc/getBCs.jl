@views function getBCs_D(Abc,A,nx,ny,nD,BCtype)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    # calculation
    for dim ∈ 1:nD
        if BCtype == 1 #"periodic"
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = A[nx ,j,dim]
                elseif i == nx+2
                    Abc[nx+2,j,dim] = A[1  ,j,dim]
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end

        elseif BCtype == 2 #"dirichlet"
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = A[1 ,j,dim]
                elseif i == nx+2
                    Abc[nx+2,j,dim] = A[nx ,j,dim]
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end 
        elseif BCtype == 3 #"reflective"
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = A[1 ,j,dim]
                elseif i == nx+2
                    Abc[nx+2,j,dim] = A[nx ,j,dim]
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end 
        elseif BCtype == 4 #"outflow"
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = 0.0
                elseif i == nx+2
                    Abc[nx+2,j,dim] = 0.0
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end 
        end
    end
    return nothing
end

@views function getBCs(A,nx,ny,nD,type)
    Abc = zeros(nx+2,ny,nD)
    if type == "periodic"
        for dim ∈ 1:nD
            for j ∈ 1:ny
                for i ∈ 1:nx+2
                    if i == 1
                        Abc[1   ,j,dim] = A[nx ,j,dim]
                    elseif i == nx+2
                        Abc[nx+2,j,dim] = A[1  ,j,dim]
                    else
                        Abc[i   ,j,dim] = A[i-1,j,dim]
                    end
                end
            end
        end
    elseif type == "dirichlet"
        for dim ∈ 1:nD
            for j ∈ 1:ny
                for i ∈ 1:nx+2
                    if i == 1
                        Abc[1   ,j,dim] = A[1   ,j,dim]
                    elseif i == nx+2
                        Abc[nx+2,j,dim] = A[nx,j,dim]
                    else
                        Abc[i   ,j,dim] = A[i-1 ,j,dim]
                    end
                end
            end
        end  
    elseif type == "reflective"
            for j ∈ 1:ny
                for i ∈ 1:nx+2
                    if i == 1
                        Abc[1   ,j,1] =  A[1   ,j,1]
                        Abc[1   ,j,2] = -A[1   ,j,2]
                        Abc[1   ,j,3] = -A[1   ,j,3]
                    elseif i == nx+2
                        Abc[nx+2,j,1] =  A[nx,j,1]
                        Abc[nx+2,j,2] = -A[nx,j,2]
                        Abc[nx+2,j,3] = -A[nx,j,3]
                    else
                        Abc[i   ,j,1] =  A[i-1 ,j,1]
                        Abc[i   ,j,2] =  A[i-1 ,j,2]
                        Abc[i   ,j,3] =  A[i-1 ,j,3]
                    end
                end
            end
    elseif type == "outflow"
        for j ∈ 1:ny
            for i ∈ 1:nx+2
                if i == 1
                    Abc[1   ,j,1] = 0.0
                    Abc[1   ,j,2] = 0.0
                    Abc[1   ,j,3] = 0.0
                elseif i == nx+2
                    Abc[nx+2,j,1] = 0.0
                    Abc[nx+2,j,2] = 0.0
                    Abc[nx+2,j,3] = 0.0
                else
                    Abc[i   ,j,1] =  A[i-1 ,j,1]
                    Abc[i   ,j,2] =  A[i-1 ,j,2]
                    Abc[i   ,j,3] =  A[i-1 ,j,3]
                end
            end
        end
    end
    return Abc
end

#=
@views function getBCs(A,nx,ny,nD)
    Abc = zeros(nx+2,ny,nD)
    for dim in 1:nD
        for j in 1:ny
            for i in 1:nx+2
                if i == 1
                    Abc[i,j,dim] = 0.0*A[2 ,j,dim]
                elseif i == nx+2
                    Abc[i,j,dim] = 0.0*A[nx-1,j,dim]
                else
                    Abc[i,j,dim] = A[i-1,j,dim]
                end
            end
        end
    end
    return Abc
end
=#