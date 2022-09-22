@views function advSolve_D(cublocks,cuthreads,h,Qx,Qy,UFS,Ubc,U,zbc,z,g,Δx,Δy,Δt,nx,ny,type)
    # x-direction, F flux vector
    # assembly of conservative variables vector and flux function vector
        @cuda blocks=cublocks threads=cuthreads getU_D(U,h,Qx,Qy,nx,ny)
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc,Ubc,z,U,nx,ny,1)
    # get fluxes x-direction
        if type=="Rus"
            @cuda blocks=cublocks threads=cuthreads fluxRus_D(UFS,Ubc,zbc,g,nx,ny,1)
        elseif type=="HLL"
            @cuda blocks=cublocks threads=cuthreads fluxHLL_D(UFS,Ubc,zbc,g,nx,ny,1)
        elseif type=="HLLC"
            @cuda blocks=cublocks threads=cuthreads fluxHLLC_D(UFS,Ubc,zbc,g,nx,ny,1)
        else 
            @error "invalid numerical flux definition, valid ones are:\n\t a) Rus  - Rusanov fluxes\n\t b) HLL  - HLL approximate Riemann solver\n\t c) HLLC - HLLC  approximate Riemann solver"
            exit(-1)
        end
    # update along x-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U,UFS,(Δt/Δx),nx,ny,1)
    synchronize()

    # y-direction, G flux vector
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc,Ubc,z,U,nx,ny,2)
    # get fluxes y-direction
    if type=="Rus"
        @cuda blocks=cublocks threads=cuthreads fluxRus_D(UFS,Ubc,zbc,g,nx,ny,2)
    elseif type=="HLL"
        @cuda blocks=cublocks threads=cuthreads fluxHLL_D(UFS,Ubc,zbc,g,nx,ny,2)
    elseif type=="HLLC"
        @cuda blocks=cublocks threads=cuthreads fluxHLLC_D(UFS,Ubc,zbc,g,nx,ny,2)
    else 
        @error "invalid numerical flux definition, valid ones are:\n\t a) Rus  - Rusanov fluxes\n\t b) HLL  - HLL approximate Riemann solver\n\t c) HLLC - HLLC  approximate Riemann solver"
        exit(-1)
    end
    # update along y-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U,UFS,(Δt/Δy),nx,ny,2)
        @cuda blocks=cublocks threads=cuthreads getQxQyh_D(h,Qx,Qy,U,g,nx,ny)
    synchronize()
end