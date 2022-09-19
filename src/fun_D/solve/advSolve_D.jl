# include dependencies & function call(s)
include("get_D.jl")
include("../bc/getBCs_D.jl")
include("../flux/fluxes_D.jl")
include("../upd/update_D.jl")
@views function advSolve_D(cublocks,cuthreads,h,Qx,Qy,UFS,Ubc,U,zbc,z,g,Δx,Δy,Δt,nx,ny,type)
    # x-direction, F flux vector
    # assembly of conservative variables vector and flux function vector
        @cuda blocks=cublocks threads=cuthreads getU_D(U,h,Qx,Qy,nx,ny)
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc,Ubc,z,U,nx,ny,1)
    # get fluxes x-direction
        @cuda blocks=cublocks threads=cuthreads fluxRus_D(UFS,Ubc,zbc,g,nx,ny,1) 
    # update along x-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U,UFS,(Δt/Δx),nx,ny,1)
    synchronize()

    # y-direction, G flux vector
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc,Ubc,z,U,nx,ny,2)
    # get fluxes y-direction
        @cuda blocks=cublocks threads=cuthreads fluxRus_D(UFS,Ubc,zbc,g,nx,ny,2) 
    # update along y-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U,UFS,(Δt/Δy),nx,ny,2)
        @cuda blocks=cublocks threads=cuthreads getQxQyh_D(h,Qx,Qy,U,g,nx,ny)
    synchronize()
end