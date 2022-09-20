@views function getBCs_D(Abc,A,nx,ny,nD,BCtype)
    # index initialization
    i = (blockIdx().x-1)*blockDim().x+threadIdx().x
    j = (blockIdx().y-1)*blockDim().y+threadIdx().y
    # calculation
    for dim ∈ 1:nD
        if BCtype == 1 #periodic
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = A[nx ,j,dim]
                elseif i == nx+2
                    Abc[nx+2,j,dim] = A[1  ,j,dim]
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end

        elseif BCtype == 2 #dirichlet
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = A[1 ,j,dim]
                elseif i == nx+2
                    Abc[nx+2,j,dim] = A[nx ,j,dim]
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end 
        elseif BCtype == 3 #reflective
            if i<=nx && j<=ny
                if i == 1
                    Abc[1   ,j,dim] = A[1 ,j,dim]
                elseif i == nx+2
                    Abc[nx+2,j,dim] = A[nx ,j,dim]
                else
                    Abc[i   ,j,dim] = A[i-1,j,dim]
                end
            end 
        elseif BCtype == 4 #outflow
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

@views function getBC_D(zbc,Ubc,z,U,nx,ny,direction)
    # index initialization
    i = (blockIdx().x-1)*blockDim().x+threadIdx().x
    j = (blockIdx().y-1)*blockDim().y+threadIdx().y
    # calculation
    if direction == 1
        for dim ∈ 1:3
            # outflow
            if i<=nx && j<=ny
                if i == 1
                    Ubc[1   ,j,dim] = 0.0
                elseif i == nx+2
                    Ubc[nx+2,j,dim] = 0.0
                else
                    Ubc[i   ,j,dim] = U[i-1,j,dim]
                end
            end 
        end
        # dirichlet
        if i<=nx && j<=ny
            if i == 1
                zbc[1   ,j] = z[1  ,j]
            elseif i == nx+2
                zbc[nx+2,j] = z[nx ,j]
            else
                zbc[i   ,j] = z[i-1,j]
            end
        end 
    elseif direction == 2
        for dim ∈ 1:3
            # outflow
            if i<=nx && j<=ny
                if j == 1
                    Ubc[i,1   ,dim] = 0.0
                elseif j == ny+2
                    Ubc[i,ny+2,dim] = 0.0
                else
                    Ubc[i,j   ,dim] = U[i,j-1,dim]
                end
            end 
        end
        # dirichlet
        if i<=nx && j<=ny
            if j == 1
                zbc[i,1   ] = z[i,1  ]
            elseif j == ny+2
                zbc[i,ny+2] = z[i,ny ]
            else
                zbc[i   ,j] = z[i,j-1]
            end
        end 
    end
    return nothing
end