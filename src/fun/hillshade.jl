@views function hillshade(z,Δx,Δy,ϕ,θ,nx,ny)
    #=
    -------------
    | a | b | c |
    -------------
    | d | e | f | e = i,j
    -------------
    | g | h | i |
    -------------
    =#
    as     = zeros(Float64,nx,ny)
    hs     = zeros(Float64,nx,ny)
    for j ∈ 2:ny-1
        for i ∈ 2:nx-1
            A = z[i-1,j-1]
            B = z[i-1,j  ]
            C = z[i-1,j+1]
            D = z[i  ,j-1]
            F = z[i  ,j+1]
            G = z[i+1,j-1]
            H = z[i+1,j  ]
            I = z[i+1,j+1]

            ∂zx = ((C+2.0*F+I)-(A+2.0*D+G))/(8.0*Δx)
            ∂zy = ((G+2.0*H+I)-(A+2.0*B+C))/(8.0*Δy)
            s   = atand(sqrt(∂zx^2+∂zy^2))
            a   = atand(∂zy,-∂zx)
            
            if a<90.0
                as[i,j] = -a+90.0
            elseif a>=90.0
                as[i,j] = 360.0-a+90.0
            end
            h = 255.0*((cosd(ϕ)*cosd(s))+(sind(ϕ)*sind(s*cosd(θ-as[i,j]))))
            if h >= 0.0
                hs[i,j] = abs(h+1.0)
            else
                hs[i,j] = 1.0
            end
        end
    end
    return hs
end