@views function updateU!(U,F,c,nx,ny,nD)
    for dim ∈ 1:nD
        for j ∈ 1:ny
            for i ∈ 1:nx
                U[i,j,dim] -= c*F[i,j,dim]
            end
        end
    end
    return nothing
end
@views function updateAdvU!(U,S,Δt,nx,ny,nD,flow_type)
    if flow_type=="coulomb"
        for dim ∈ 1:nD
            for j ∈ 1:ny
                for i ∈ 1:nx
                    U[i,j,dim]+=Δt*S[i,j,dim]
                end
            end
        end
    end
    if flow_type=="newtonian" || flow_type=="plastic"
        for dim ∈ 1:3
            for j ∈ 1:ny
                for i ∈ 1:nx
                    if dim == 1
                        U[i,j,1]+=Δt*S[i,j,1]
                    else
                        U[i,j,dim]=U[i,j,dim]/(1.0+Δt*S[i,j,dim])
                    end
                end
            end
        end
    end
    return nothing
end