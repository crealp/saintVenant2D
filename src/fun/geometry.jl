@views function oneDkern(z,xc,yc,i,j,nx,ny,r)
    val = z[i,j]
    xc0 = xc[i]
    yc0 = yc[j]
    c   = 0
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = sqrt((xc[i]-xc0)^2+(yc[j]-yc0)^2)
            if d<=r
                val += z[i,j]*(r/(d+r))
                c   += 1
            end
        end
    end
    return val/c
    #=
    σ   = r
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = sqrt((xc[i]-xc0)^2+(yc[j]-yc0)^2)
            if d<=r
                val += z[i,j]*1.0/(sqrt(2.0*pi*σ))*exp(-d^2/(2.0*σ^2))
            end
        end
    end
    return val
    =#
end

@views function geometry(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    z      = 0.0.*exp.(((xc.-lx/2)./(lx./2)).^2 .+((yc.-ly/2)./(ly./2))'.^2)
    δz     = 0.0.*1.5.*(2.0.*rand(nx,ny).-1.0)
    for r ∈ 1:5
        for i ∈ 1:nx
            for j ∈ 1:ny
                if i>1 && i<nx && j>1 && j<ny
                    δz[i,j] = (δz[i-1,j-1]+δz[i,j-1]+δz[i+1,j-1]+δz[i-1,j]+δz[i,j]+δz[i+1,j]+δz[i-1,j+1]+δz[i,j+1]+δz[i+1,j+1])/9.0
                elseif i==1
                    δz[i,j] = δz[i+1,j  ]
                elseif i==nx
                    δz[i,j] = δz[i-1,j  ]
                elseif j==1
                    δz[i,j] = δz[i  ,j+1]
                elseif j==ny
                    δz[i,j] = δz[i  ,j-1]    
                end

            end
        end
    end
    z = z.+δz
    #z = (tanh.((yc.-0.75*ly)./15).*ones(Float64,1,nx))'
    # set initial fluid height
    hi = 1.0/2.0 
    #hi = 2.0
    h0 = hi.-(z)
    h  = h0.+exp.(.-((xc.-lx/3)./(lx./8)).^2 .-((yc.-lx/2)./(lx./8))'.^2)
    h  = h0.+exp.(.-((xc.-lx/2)./(lx./8)).^2 .-((yc.-lx/8)./(lx./4))'.^2)
    
    return(h,z,xc,yc,Δx,Δy)
end

@views function gaussian_floor(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    z      = 0.0.*exp.(.-((xc.-lx/2)./(lx./8)).^(2) .-((yc.-lx/2)./(lx./8))'.^(2))
    Δz     = 0.5
    δz     = 1.0.*Δz.*(2.0.*rand(nx,ny).-1.0)
    for r ∈ 1:2
        for i ∈ 1:nx
            for j ∈ 1:ny
                δz[i,j] = oneDkern(δz,xc,yc,i,j,nx,ny,2.25*Δx)
            end
        end
    end
    δz = (Δz/maximum(abs.(δz))).*δz


    xc0    = 0.5*lx
    yc0    = 0.5*ly
    R      = lx/20
    hbump  = 2.5
    H      = hbump/exp(-1.0/R^2)
    for i ∈ 1:nx
        for j ∈ 1:ny
            r = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if r<=R
                z[i,j]+=H*exp(-1.0/(R^2-((xc[i]-xc0)^2+(yc[j]-yc0)^2)  ))

            end
        end
    end

    z = z.+δz
    # set initial fluid height
    hi = 1.5
    h0 = hi.-z
    h  = h0.*ones(Float64,nx,ny)
    h  = h0.+exp.(.-((xc.-lx/2)./(lx./16)).^(2) .-((yc.-lx/4)./(lx./16))'.^(2))
    
    return(h,z,xc,yc,Δx,Δy)
end

@views function hyperbolic_floor(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    Δz     = 0.5
    δz     = 1.0.*Δz.*(2.0.*rand(nx,ny).-1.0)
    nr     = 2
    n      = nr+nx*ny
    println("[=> generating randomly distributed floor...")
    p      = Progress(n, dt=0.1, barglyphs=BarGlyphs("[=> ]"), barlen=10, color=:green, showspeed=true)
    for r ∈ 1:nr
        for i ∈ 1:nx
            for j ∈ 1:ny
                δz[i,j] = oneDkern(δz,xc,yc,i,j,nx,ny,0.5)
                next!(p)
            end
        end
    end
    δz = (Δz/maximum(abs.(δz))).*δz
    z  = (tanh.((yc.-0.75*ly)./0.5).*ones(Float64,1,nx))'
    z  = z.+δz
    # set initial fluid height
    hi = 1.0 
    h  = (hi.-z).+exp.(.-((xc.-lx/2)./(min(lx,ly)./16)).^(2) .-((yc.-ly/2)./(min(lx,ly)./16))'.^(2))
    
    return(h,z,xc,yc,Δx,Δy)
end
@views function bowl_floor(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    z      = exp.(((xc.-lx/2)./(lx./3)).^2 .+((yc.-ly/2)./(ly./3))'.^2)
    Δz     = 0.5
    δz     = 1.0.*Δz.*(2.0.*rand(nx,ny).-1.0)
    nr     = 2
    n      = nr+nx*ny
    println("[=> generating randomly distributed floor...")
    p      = Progress(n, dt=0.1, barglyphs=BarGlyphs("[=> ]"), barlen=10, color=:green, showspeed=true)
    for r ∈ 1:nr
        temp = zeros(Float64,nx,ny)
        for j ∈ 1:ny
            for i ∈ 1:nx
                temp[i,j] = oneDkern(δz,xc,yc,i,j,nx,ny,min(lx,ly)/10.0)
                next!(p)
            end
        end
        δz = temp
        δz = (Δz/maximum(abs.(δz))).*δz
    end
    z = z.+δz
    #z = (tanh.((yc.-0.75*ly)./15).*ones(Float64,1,nx))'
    # set initial fluid height
    hi = 1.0/2.0
    h  = zeros(Float64,nx,ny)
    xc0 = lx/2
    yc0 = 0.2*ly
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if d<=lx/10
                h[i,j]=hi-δz[i,j]
            end
        end
    end
    xc0 = lx/2
    yc0 = 0.8*ly
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if d<=lx/10
                h[i,j]=hi-δz[i,j]
            end
        end
    end
    return(h,z,xc,yc,Δx,Δy)
end
@views function bowl_floor_noΔz(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    z      = exp.(((xc.-lx/2)./(lx./3.0)).^2 .+((yc.-ly/2)./(ly./3.0))'.^2)
    # set initial fluid height
    hi = 1.0/2.0
    h  = zeros(Float64,nx,ny)
    xc0 = lx/2
    yc0 = 0.2*ly
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if d<=lx/10
                h[i,j]=hi
            end
        end
    end
    xc0 = lx/2
    yc0 = 0.8*ly
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if d<=lx/10
                h[i,j]=hi
            end
        end
    end
#=
    h  = zeros(Float64,nx,ny)
    xc0 = lx/2
    yc0 = lx/2
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if d>=lx/1.5 && d<=lx/1.2
                h[i,j]=hi
            end
        end
    end
=#
    return(h,z,xc,yc,Δx,Δy)
end
@views function bowl_floor_rough(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    z      = exp.(((xc.-lx/2)./(lx./2)).^2 .+((yc.-ly/2)./(ly./2))'.^2)
    Δz     = 3.5
    δz     = 1.0.*Δz.*(2.0.*rand(nx,ny).-1.0)
    nr     = 2
    n      = nr+nx*ny
    println("[=> generating randomly distributed floor...")
    p      = Progress(n, dt=0.1, barglyphs=BarGlyphs("[=> ]"), barlen=10, color=:green, showspeed=true)
    for r ∈ 1:nr
        temp = zeros(Float64,nx,ny)
        for j ∈ 1:ny
            for i ∈ 1:nx
                temp[i,j] = oneDkern(δz,xc,yc,i,j,nx,ny,2.0)
                next!(p)
            end
        end
        δz = temp
        δz = (Δz/maximum(abs.(δz))).*δz
    end
    z = z.+δz
    #z = (tanh.((yc.-0.75*ly)./15).*ones(Float64,1,nx))'
    # set initial fluid height
    hi = 1.0e-3
    h  = hi.*ones(Float64,nx,ny)
    return(h,z,xc,yc,Δx,Δy)
end
@views function incline(lx,ly,nx,ny)
    # number of points
    Δx,Δy  = lx/nx,ly/ny 
    x,y    = 0.0:Δx:lx,0.0:Δy:ly
    # calculate midpoint values of x in each control vlume
    xc,yc  = 0.5.*(x[1:nx]+x[2:nx+1]),0.5.*(y[1:ny]+y[2:ny+1])
    # set initial bed topography
    a      = -1.0
    zmax   = 0.5*(abs(a)*lx)
    z      = exp.(((xc.-lx/2)./(lx./2)).^2 .+((yc.-ly/2)./(ly./2))'.^2)
    z      = (a.*xc.+zmax).*ones(Float64,ny)'
    for j ∈ 1:ny
        for i ∈ 1:nx
            if z[i,j]<=-1e-10
                z[i,j]=-1e-10
            end
        end
    end
    xc0    = lx/2.5
    yc0    = 0.5*ly
    R      = lx/20
    hbump  = 0.91
    H      = hbump/exp(-1.0/R^2)
    for i ∈ 1:nx
        for j ∈ 1:ny
            r = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if r<=R
                z[i,j]+=H*exp(-1.0/(R^2-((xc[i]-xc0)^2+(yc[j]-yc0)^2)  ))

            end
        end
    end
    # set initial fluid height
    hi     = 1.0/2.0
    h      = zeros(Float64,nx,ny)
    xc0    = lx/5
    yc0    = 0.5*ly
    for i ∈ 1:nx
        for j ∈ 1:ny
            d = (xc[i]-xc0)^2+(yc[j]-yc0)^2
            if d<=lx/10
                h[i,j]=max((a.*xc0.+zmax+hi)-z[i,j],0.0)

            end
        end
    end
    return(h,z,xc,yc,Δx,Δy)
end