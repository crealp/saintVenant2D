# include dependencies & function call(s)
using Plots, LaTeXStrings, Base.Threads, ProgressMeter, CSV, DataFrames, CUDA

@views function hillshade_D!(as,hs,z,Δx,Δy,ϕ,θ,nx,ny)
    # index initialization
    i  = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j  = (blockIdx().y-1) * blockDim().y + threadIdx().y
    # pre-computation
    Δx = 1.0/(8.0*Δx)
    Δy = 1.0/(8.0*Δy)
    if i>1 && i<=(nx-1) && j>1 && j<=(ny-1)
        A   = z[i-1,j-1]
        B   = z[i-1,j  ]
        C   = z[i-1,j+1]
        D   = z[i  ,j-1]
        F   = z[i  ,j+1]
        G   = z[i+1,j-1]
        H   = z[i+1,j  ]
        I   = z[i+1,j+1]

        ∂zx = ((C+2.0*F+I)-(A+2.0*D+G))*Δx
        ∂zy = ((G+2.0*H+I)-(A+2.0*B+C))*Δy
        #= TO BE TESTED
        ∂zx = ((z[i-1,j+1]+2.0*z[i  ,j+1]+z[i+1,j+1])-(z[i-1,j-1]+2.0*z[i  ,j-1]+z[i+1,j-1]))*Δx
        ∂zy = ((z[i+1,j-1]+2.0*z[i+1,j  ]+z[i+1,j+1])-(z[i-1,j-1]+2.0*z[i-1,j  ]+z[i-1,j+1]))*Δy
        =#
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
    return nothing
end
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
    path = "./docs/py_vs_ju/dat/parameters.csv"
    # get infos in headers
    D = Array(CSV.read(path,DataFrame,header=true; delim=",")[1,:])
    nx= Int64(D[1])
    ny= Int64(D[2])
    dx= Float64(D[3])
    dy= Float64(D[4])
    # get infos in headers
    path = "./docs/py_vs_ju/dat/x.csv"
    xc= Array(CSV.read(path,DataFrame,header=1; delim=",")[:,1])
    path = "./docs/py_vs_ju/dat/y.csv"
    yc= Array(CSV.read(path,DataFrame,header=1; delim=",")[:,1])
    # get infos in headers
    path = "./docs/py_vs_ju/dat/zhs.csv"
    z = reshape(Array(CSV.read(path,DataFrame,header=1; delim=",")[:,1]),nx,ny)


    #z = rand(Float64,21000,21000)
    #nx,ny = size(z)

    # define grid & block sizes for kernel initialization
    BLOCKx    = 8
    BLOCKy    = 8
    GRIDx     = ceil(Int,nx/BLOCKx)
    GRIDy     = ceil(Int,ny/BLOCKy)
    cuthreads = (BLOCKx, BLOCKy, 1)
    cublocks  = (GRIDx,  GRIDy,  1)

    @info "GPU kernel:" cuthreads,cublocks nx,ny
    # allocate memory on GPU, i.e., variable on device specified by <nameOfVariable>_D
    z_D  = CUDA.zeros(Float64,nx,ny)
    as_D = CUDA.zeros(Float64,nx,ny)
    hs_D = CUDA.zeros(Float64,nx,ny)
    # copy from CPU to GPU
    copyto!(z_D,z)
    for k in 1:5
        CUDA.@time @cuda blocks=cublocks threads=cuthreads hillshade_D!(as_D,hs_D,z_D,dx,dy,45.0,315.0,nx,ny)
        synchronize()
    end
    println()
    for k in 1:5
        @time hs=hillshade(z,dx,dy,45.0,315.0,nx,ny)
    end
    hs = Array(hs_D)
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        hillshade_plot(xc,yc,hs,45.0,315.0,0.75)
    savefig("./docs/py_vs_ju/julia_hillshade.png")
    

end
main()