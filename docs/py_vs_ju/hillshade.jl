# include dependencies & function call(s)
using Plots, LaTeXStrings, Base.Threads, ProgressMeter, CSV, DataFrames


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

@views function main()
    path = "./scripts/py_vs_ju/dat/parameters.csv"
    # get infos in headers
    D = CSV.read(path,DataFrame,header=true; delim=",")
    D = Array(D[1,:])
    nx= Int64(D[1])
    ny= Int64(D[2])
    dx= Float64(D[3])
    dy= Float64(D[4])
    path = "./scripts/py_vs_ju/dat/x.csv"
    # get infos in headers
    D = CSV.read(path,DataFrame,header=1; delim=",")
    xc= Array(D[:,1])
    path = "./scripts/py_vs_ju/dat/y.csv"
    # get infos in headers
    D = CSV.read(path,DataFrame,header=1; delim=",")
    yc= Array(D[:,1])
    path = "./scripts/py_vs_ju/dat/zhs.csv"
    # get infos in headers
    D = CSV.read(path,DataFrame,header=1; delim=",") 
    z = reshape(Array(D[:,1]),nx,ny)
    @time hs=hillshade(z,dx,dy,45.0,315.0,nx,ny)
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        hillshade_plot(xc,yc,hs,45.0,315.0,0.75)
    savefig("./scripts/py_vs_ju/julia_hillshade.png")
end
main()