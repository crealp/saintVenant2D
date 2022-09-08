include("../upd/update.jl")
@views function τ_coulomb(S,h,Qx,Qy,z,g,nx,ny,Δx,Δy)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y

    ρs = 2.7e3          # solid density
    ϕb = 35.0*pi/180    # internal friction angle
    μ0 = tan(ϕb)        # static friction coefficient
    μw = tan(0.5*ϕb)    # dynamic friction coefficient
    W  = 1.0e6          # velocity threshold
    if i<=nx && j<=ny
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
    return nothing
end

@views function souSolve_D(cublocks,cuthreads,h,Qx,Qy,z,U,g,Δx,Δy,t,Δt,nx,ny,flow_type,pcpt_onoff)
    S  = CUDA.zeros(Float64,nx,ny,3)
    @cuda blocks=cublocks threads=cuthreads τ_coulomb(S,h,Qx,Qy,z,g,nx,ny,Δx,Δy)
    synchronize()
    # assembly of conservative variables vector and flux function vector
    @cuda blocks=cublocks threads=cuthreads getU_D(U,h,Qx,Qy,nx,ny)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads updateAdvU_D(U,S,Δt,nx,ny)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads getQxQyh_D(h,Qx,Qy,U,g,nx,ny)
    synchronize()
    return nothing
end