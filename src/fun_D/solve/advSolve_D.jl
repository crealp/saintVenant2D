# include dependencies & function call(s)
include("get_D.jl")
include("../bc/getBCs_D.jl")
include("../flux/fluxes_D.jl")
include("../upd/update_D.jl")
@views function advSolve_D(cublocks,cuthreads,h_D,Qx_D,Qy_D,UFS,Ubc_D,U_D,zbc_D,z_D,g,Δx,Δy,Δt,nx,ny,type)
    # x-direction, F flux vector
    # assembly of conservative variables vector and flux function vector
        @cuda blocks=cublocks threads=cuthreads getU_D(U_D,h_D,Qx_D,Qy_D,nx,ny)
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc_D,Ubc_D,z_D,U_D,nx,ny,1)
    # get fluxes x-direction
        @cuda blocks=cublocks threads=cuthreads fluxRus2_D(UFS,Ubc_D,zbc_D,g,nx,ny,1) 
    # update along x-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U_D,UFS,(Δt/Δx),nx,ny,1)
    synchronize()

    # y-direction, G flux vector
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc_D,Ubc_D,z_D,U_D,nx,ny,2)
    # get fluxes y-direction
        @cuda blocks=cublocks threads=cuthreads fluxRus2_D(UFS,Ubc_D,zbc_D,g,nx,ny,2) 
    # update along y-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U_D,UFS,(Δt/Δy),nx,ny,2)
        @cuda blocks=cublocks threads=cuthreads getQxQyh_D(h_D,Qx_D,Qy_D,U_D,g,nx,ny)
    synchronize()
end