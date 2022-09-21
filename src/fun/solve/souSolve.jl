@views function τ_coulomb!(S,h,Qx,Qy,z,g,nx,ny,Δx,Δy)
    ρs = 2.7e3          # solid density
    ϕb = 35.0*pi/180    # internal friction angle
    μ0 = tan(ϕb)        # static friction coefficient
    μw = tan(0.5*ϕb)    # dynamic friction coefficient
    W  = 1.0e6          # velocity threshold
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j]>0.0
                u  = Qx[i,j]/(h[i,j])   # x-component velocity
                v  = Qy[i,j]/(h[i,j])   # y-component velocity
                w  = sqrt(u^2+v^2)      # magnitude L2 of the velocity
                if i==1
                    αx = atan((z[i+1,j]-z[nx ,j])/(2.0*Δx))                            
                elseif i==nx
                    αx = atan((z[1  ,j]-z[i-1,j])/(2.0*Δx))
                else 
                    αx = atan((z[i+1,j]-z[i-1,j])/(2.0*Δx))                            
                end
                if j==1
                    αy = atan((z[i,j+1]-z[i,ny ])/(2.0*Δy))                            
                elseif j==ny
                    αy = atan((z[j,1  ]-z[i,j-1])/(2.0*Δy))
                else 
                    αy = atan((z[i,j+1]-z[i,j-1])/(2.0*Δy))                            
                end
                if w>0.0
                    μ  = (μ0-μw)/(1.0+w/W)+μw   # velocity-dependent friction model, see yamada etal, 2018
                    τ  = ρs*g*h[i,j]*μ          # basal frictional/shear resistance law, see 
                    τx = τ*cos(αx)*(u/w)        # x-component basal shear
                    τy = τ*cos(αy)*(v/w)        # y-component basal shear
                else 
                    τx = 0.0
                    τy = 0.0
                end
                S[i,j,1] = 0.0
                S[i,j,2] = -τx/ρs
                S[i,j,3] = -τy/ρs                   
            end
        end
    end
    return S
end
@views function τ_newtonian(h,Qx,Qy,g,nx,ny)
    n  = 0.00025
    ρw = 1.0e3
    S  = zeros(Float64,nx,ny,3)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j]>0.0
                u  = Qx[i,j]/(h[i,j])
                v  = Qy[i,j]/(h[i,j])
                w  = sqrt(u^2+v^2)
                if w>0.0
                        Cf = n^2*(h[i,j])^(-4/3)
                        τ  = g*Cf*w
                else 
                        τ  = 0.0
                end
                S[i,j,1] = 0.0
                S[i,j,2] = τ
                S[i,j,3] = τ
            end
        end
    end
    return S
end
@views function τ_plastic(h,Qx,Qy,g,nx,ny)
    ρs = 2.7e3
    ϕb = 15.0*pi/180
    μ  = tan(ϕb)
    η  = 1.0e3
    m  = 3.0/2.0
    S  = zeros(Float64,nx,ny,3)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j]>0.0
                u  = Qx[i,j]/(h[i,j])
                v  = Qy[i,j]/(h[i,j])
                w  = sqrt(u^2+v^2)
                u  = Qx[i,j]/(h[i,j])
                v  = Qy[i,j]/(h[i,j])
                w  = sqrt(u^2+v^2)

                τf = ρs*g*h[i,j]*tan(ϕb)
                

                τηx = ((2.0*m+1.0)/m)^m*η*(abs(u)/h[i,j])^m
                τηy = ((2.0*m+1.0)/m)^m*η*(abs(v)/h[i,j])^m


                S[i,j,1] = 0.0
                S[i,j,2] = (τf+τηx)/ρs
                S[i,j,3] = (τf+τηy)/ρs
            end
        end
    end
    return S
end
@views function precip!(S,ϵp,t,nx,ny)
    A = (2.0*rand(Float64)-1.0)
    f = 1.0
    r = A*sin(2*pi*5.0*t)
    p = 0.0

    for j ∈ 1:ny
        for i ∈ 1:nx
            S[i,j,1] = ϵp
        end
    end
    return nothing
end
@views function souSolve(h,Qx,Qy,z,U,g,Δx,Δy,t,Δt,nx,ny,flow_type,pcpt_onoff)
    S  = zeros(Float64,nx,ny,3)
    if flow_type=="coulomb"
        τ_coulomb!(S,h,Qx,Qy,z,g,nx,ny,Δx,Δy)
    elseif flow_type=="newtonian"
        S = τ_newtonian(h,Qx,Qy,g,nx,ny)
    elseif flow_type=="plastic"
        S = τ_plastic(h,Qx,Qy,g,nx,ny)
    end
    # add precipitation if pcpt_onoff==true
    if pcpt_onoff==true
        precip!(S,8.0e-6,t,nx,ny)
    end
    # assembly of conservative variables vector and flux function vector
    getU!(U,h,Qx,Qy,nx,ny)
    updateAdvU!(U,S,Δt,nx,ny,3,flow_type)
    return copy(U[:,:,1]),copy(U[:,:,2]),copy(U[:,:,3])
end