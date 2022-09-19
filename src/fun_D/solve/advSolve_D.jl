# include dependencies & function call(s)
include("get_D.jl")
include("../bc/getBCs_D.jl")
include("../flux/fluxes_D.jl")
include("../upd/update_D.jl")
@views function advSolve_D(cublocks,cuthreads,h_D,Qx_D,Qy_D,UFS,Ubc_D,U_D,zbc_D,z_D,g,Δx,Δy,Δt,nx,ny,type)
    # assembly of conservative variables vector and flux function vector
        @cuda blocks=cublocks threads=cuthreads getU_D(U_D,h_D,Qx_D,Qy_D,nx,ny)
        @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
    # ghost cells
        @cuda blocks=cublocks threads=cuthreads getBC_D(zbc_D,Ubc_D,z_D,U_D,nx,ny,1)
    # get fluxes
        @cuda blocks=cublocks threads=cuthreads fluxRus_D(UFS,Ubc_D,zbc_D,g,nx,ny,1) 
    # update along x-direction
        @cuda blocks=cublocks threads=cuthreads updateU_D(U_D,UFS,(Δt/Δx),nx,ny,1)
    synchronize()

    # ghost cells
    U_D   = permutedims(U_D  ,(2,1,3))
    z_D   = permutedims(z_D  ,(2,1)  )
    Ubc_D = permutedims(Ubc_D,(2,1,3))
    zbc_D = permutedims(zbc_D,(2,1)  )
    @cuda blocks=cublocks threads=cuthreads setUFS_D(UFS,nx+1,ny+1)
    synchronize()
    UFS   = permutedims(UFS ,(2,1,3,4))
    @cuda blocks=cublocks threads=cuthreads getBC_D(zbc_D,Ubc_D,z_D,U_D,ny,nx,1)
    synchronize()

    @cuda blocks=cublocks threads=cuthreads fluxRus_D(UFS,Ubc_D,zbc_D,g,ny,nx,2) 
    synchronize()

    UFS = permutedims(UFS ,(2,1,3,4))
    U_D = permutedims(U_D ,(2,1,3))
    z_D = permutedims(z_D ,(2,1)  )
    @cuda blocks=cublocks threads=cuthreads updateU_D(U_D,UFS,(Δt/Δy),nx,ny,2)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads getQxQyh_D(h_D,Qx_D,Qy_D,U_D,g,nx,ny)
    synchronize()
end