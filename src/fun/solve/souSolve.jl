@views function updateAdvU!(U,S,Δt,nx,ny,nD,flow_type)
    if flow_type=="coulomb"
        for dim in 1:nD
            for j in 1:ny
                for i in 1:nx
                    U[i,j,dim]+=Δt*S[i,j,dim]
                end
            end
        end
    end
    if flow_type=="newtonian"||flow_type=="plastic"
        for dim in 1:3
            for j in 1:ny
                for i in 1:nx
                U[i,j,dim]=U[i,j,dim]/(1.0-Δt*S[i,j,dim])
                end
            end
        end
    end
end
@views function precip(ϵp,t)
    A = (2.0*rand(Float64)-1.0)
    f = 1.0
    r = A*sin(2*pi*5.0*t)
    p = 0.0
    if r>0.0 && t<3600.0
        p=r*ϵp
    end
    return p
end
@views function souSolve(h,Qx,Qy,z,U,F,G,g,Δx,Δy,t,Δt,nx,ny,flow_type,pcpt_onoff)
    # assembly of conservative variables vector and flux function vector
    n  = 0.25
    m  = 1.0
    ρw = 1.0e3
    ρs = 2.7e3
    ϕb = 15.0*pi/180
    μ  = tan(ϕb)
    S  = zeros(Float64,nx,ny,3)
    for j ∈ 1:ny
        for i ∈ 1:nx
            if h[i,j]>0.0
                u  = Qx[i,j]/(h[i,j])
                v  = Qy[i,j]/(h[i,j])
                w  = sqrt(u^2+v^2)
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
                if flow_type=="coulomb"
                    if w>0.0
                        τ  = ρs*g*h[i,j]*tan(ϕb)
                        τx = τ*cos(αx)*(u/w)
                        τy = τ*cos(αy)*(v/w)
                    else 
                        τx = 0.0
                        τy = 0.0
                    end
                    S[i,j,1] = 0.0
                    S[i,j,2] = -τx/ρs
                    S[i,j,3] = -τy/ρs                   
                elseif flow_type=="plastic"
                    τ0 = 1.25e3
                    τx = 0.0
                    τy = 0.0
                    if u!=0.0
                        τx = τ0*u/w
                    end
                    if v!=0.0
                        τy = τ0*v/w
                    end
                    S[i,j,1] = 0.0
                    S[i,j,2] = -τx
                    S[i,j,3] = -τy
                elseif flow_type=="newtonian"                
                    if w>0.0
                        Cf = n^2*(h[i,j])^(-4/3)
                        τ  = g*h[i,j]*Cf*w
                    else 
                        τ  = 0.0
                    end
                    S[i,j,1] = 0.0
                    S[i,j,2] = -τ/ρw
                    S[i,j,3] = -τ/ρw
                end
            end
            if pcpt_onoff==true
                p        = precip(1.0e-3/3600.0,t)
                S[i,j,1] = p
            end
        end
    end 
    getU!(U,h,Qx,Qy,nx,ny)
    updateAdvU!(U,S,Δt,nx,ny,3,flow_type)
    return copy(U[:,:,1]),copy(U[:,:,2]),copy(U[:,:,3])
end