@views function hillshade(xc,yc,z,Δx,Δy,ϕ,θ,nx,ny)
    #=
    -------------
    | a | b | c |
    -------------
    | d | e | f | e = i,j
    -------------
    | g | h | i |
    -------------
    =#
    ∇z     = zeros(Float64,nx,ny,3)
    as     = zeros(Float64,nx,ny)
    hs     = zeros(Float64,nx,ny)
    for j ∈ 2:ny-1
        for i ∈ 2:nx-1
            A = z[i-1,j-1]
            B = z[i-1,j  ]
            C = z[i-1,j+1]
            D = z[i  ,j-1]
            E = z[i  ,j  ]
            F = z[i  ,j+1]
            G = z[i+1,j-1]
            H = z[i+1,j  ]
            I = z[i+1,j+1]

            ∇z[i,j,1] = ((C+2.0*F+I)-(A+2.0*D+G))/(8.0*Δx)
            ∇z[i,j,2] = ((G+2.0*H+I)-(A+2.0*B+C))/(8.0*Δy)
            ∇z[i,j,3] = atand(sqrt(∇z[i,j,1]^2+∇z[i,j,2]^2))
            as[i,j]   = atand(∇z[i,j,2],-∇z[i,j,1])
            
            if as[i,j]<90.0
                as[i,j] = -as[i,j]+90
            elseif as[i,j]>=90.0
                as[i,j] = 360.0-as[i,j]+90
            end
            hs[i,j]   = 255.0*((cosd(ϕ)*cosd(∇z[i,j,3]))+(sind(ϕ)*sind(∇z[i,j,3]*cosd(θ-as[i,j]))))
            if hs[i,j] >= 0.0
                hs[i,j] = abs(hs[i,j]+1.0)
            else
                hs[i,j] = 1.0
            end
        end
    end
    return hs
end