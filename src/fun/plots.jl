@views function z_plot(xc,yc,z)
    heatmap(xc, yc, z',
        c=:terrain,
        clims=(-1.0,10.0),
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        aspect_ratio=true,
        show=true,
        title=L"z(x,y)",
    )
end
@views function z_plot!(xc,yc,z)
    heatmap!(xc, yc, z',
        c=:roma,
        clims=(0.0,10.0),
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        aspect_ratio=true,
        show=true,
        title=L"z(x,y)",
    )
end
@views function h_plot(xc,yc,h,hmax,nx,ny,t,type)
    mask = ones(Float64,nx,ny)
        for j in 1:ny
            for i in 1:nx
                if h[i,j]<=1e-3
                    mask[i,j] = NaN
                end
            end
        end
    if type == "coulomb"    
    heatmap(xc, yc, (h.*mask)',
        c=cgrad(:turbid, rev=false),
        #c=cgrad(:coffee, rev=true),
        clims=(0.0,hmax),
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        aspect_ratio=true,
        show=true,
        title=L"h(x,y)"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )
    elseif type == "plastic"    
    heatmap(xc, yc, (h.*mask)',
        c=cgrad(:turbid, rev=false),
        #c=cgrad(:coffee, rev=true),
        clims=(0.0,hmax),
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        aspect_ratio=true,
        show=true,
        title=L"h(x,y)"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )
    elseif type == "newtonian"
    heatmap(xc, yc, (h.*mask)',
        c=cgrad(:Blues, rev=false),
        clims=(0.0,hmax),
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        aspect_ratio=true,
        show=true,
        title=L"h(x,y)"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )
    else 
    heatmap(xc, yc, (h.*mask)',
        c=cgrad(:viridis, rev=false),
        clims=(0.0,hmax),
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        aspect_ratio=true,
        show=true,
        title=L"h(x,y)"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )    
    end
end
@views function wave_plot(xc,yc,h,z,η0,ηmax,nx,ny,t)
    mask = ones(Float64,nx,ny)
        for j in 1:ny
            for i in 1:nx
                if h[i,j]<=1e-2
                    mask[i,j] = NaN
                end
            end
        end
    Δη = (h.+z).-η0
    heatmap(xc, yc, (Δη.*mask)',
        colorbar=true,
        #c=cgrad(:vik,rev=false),
        c=cgrad(:RdBu,rev=true),
        clims=(-ηmax,ηmax),
        aspect_ratio=true,
        show=true,
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        title=L"\Delta\eta"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )
end
@views function free_surface_plot(xc,yc,h,z,η0,Δηmax,nx,ny,t)
    mask = ones(Float64,nx,ny)
        for j in 1:ny
            for i in 1:nx
                if h[i,j]<=1e-2
                    mask[i,j] = NaN
                end
            end
        end
    η = (h.+z)
    heatmap(xc, yc, (η.*mask)',
        colorbar=true,
        #c=cgrad(:vik,rev=false),
        c=cgrad(:viridis,rev=false),
        clims=(η0-Δηmax,η0+Δηmax),
        aspect_ratio=true,
        show=true,
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        title=L"\eta(x,y)"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )
end
@views function discharge_plot(xc,yc,h,Qx,Qy,z,vmax,nx,ny,t)
    mask = ones(Float64,nx,ny)
        for j in 1:ny
            for i in 1:nx
                if h[i,j]<=1e-3
                    mask[i,j] = NaN
                end
            end
        end
        u = Qx./(h.+1e-10)
        v = Qy./(h.+1e-10)
        V = sqrt.(u.^(2).+v.^(2))               
        heatmap(xc, yc, (v.*mask)',
            colorbar=true,
            #c=cgrad(:vik,rev=false),
            c=cgrad(:RdBu,rev=true),
            clims=(-vmax,vmax),
            aspect_ratio=true,
            show=true,
            xlabel=L"x-"*"direction [m]",
            ylabel=L"y-"*"direction [m]",
            dpi=600,
            title=L"v"*" at "*L"t = "*string(round(t,digits=2))*" [s]",
        )
end
@views function profile_plot(xc,yc,h,z,zmin,η0,nx,ny,t)
    lx,ly = size(h)
    η = (h.+z)
    id = ceil(Int64,lx/2)
    plot(xc, η[:,id],
        aspect_ratio=true,
        show=false,
    )
    plot!(yc, z[id,:],
        aspect_ratio=true,
        show=true,
        ylim=(zmin,η0),
        xlabel=L"y-"*"direction [m]",
        ylabel=L"z-"*"direction [m]",
        dpi=600,
        title=L"z(y),\eta(y)"*" at "*L"t="*string(round(t,digits=2))*" [s]",
    )
end
@views function hillshade_plot(xc,yc,hs,ϕ,θ,α)
    heatmap(xc, yc, hs',
        alpha=α,
        c=:grayC,
        xlabel=L"x-"*"direction [m]",
        ylabel=L"y-"*"direction [m]",
        dpi=600,
        clims=(0.0,255.0),
        aspect_ratio=true,
        show=true,
        title="hillshade, ("*string(round(ϕ,digits=1))*", "*string(round(θ,digits=1))*")",
    )
    return hs
end