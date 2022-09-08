# include dependencies & function call(s)
include("get.jl")
include("../bc/getBCs.jl")
include("../flux/fluxes.jl")
include("../upd/update.jl")
@views function advSolve_D(cublocks,cuthreads,h_D,Qx_D,Qy_D,UFS,Ubc_D,U_D,zbc_D,z_D,g,Δx,Δy,Δt,nx,ny,type)
    # assembly of conservative variables vector and flux function vector
    @cuda blocks=cublocks threads=cuthreads getU_D(U_D,h_D,Qx_D,Qy_D,nx,ny)
    synchronize()
    
    # ghost cells
    @cuda blocks=cublocks threads=cuthreads getBCs_D(Ubc_D,U_D,nx,ny,3,4)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads getBCs_D(zbc_D,z_D,nx,ny,1,4)
    synchronize()
    UFS = CUDA.zeros(Float64,nx+1,ny,3,7)
    # (:,:,:,1) UL
    # (:,:,:,2) UR
    # (:,:,:,3) FL
    # (:,:,:,4) FR
    # (:,:,:,5) SL
    # (:,:,:,6) SR
    # (:,:,:,7) F
    @cuda blocks=cublocks threads=cuthreads wellBal_D(UFS,Ubc_D,zbc_D,g,nx,ny,1)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads Rus_D(UFS,g,nx,ny,1)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads updateU_D(U_D,UFS,(Δt/Δx),nx,ny,1)
    synchronize()

    # ghost cells
    U_D   = permutedims(U_D  ,(2,1,3))
    z_D   = permutedims(z_D  ,(2,1)  )
    Ubc_D = permutedims(Ubc_D,(2,1,3))
    zbc_D = permutedims(zbc_D,(2,1)  )
    @cuda blocks=cublocks threads=cuthreads getBCs_D(Ubc_D,U_D,ny,nx,3,4)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads getBCs_D(zbc_D,z_D,ny,nx,1,4)
    synchronize()
    UFS = CUDA.zeros(Float64,ny+1,nx,3,7)
    # (:,:,:,1) UL
    # (:,:,:,2) UR
    # (:,:,:,3) FL
    # (:,:,:,4) FR
    # (:,:,:,5) SL
    # (:,:,:,6) SR
    @cuda blocks=cublocks threads=cuthreads wellBal_D(UFS,Ubc_D,zbc_D,g,ny,nx,2)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads Rus_D(UFS,g,ny,nx,2)
    synchronize()
    UFS = permutedims(UFS ,(2,1,3,4))
    U_D = permutedims(U_D ,(2,1,3))
    z_D = permutedims(z_D ,(2,1)  )
    @cuda blocks=cublocks threads=cuthreads updateU_D(U_D,UFS,(Δt/Δy),nx,ny,2)
    synchronize()
    @cuda blocks=cublocks threads=cuthreads getQxQyh_D(h_D,Qx_D,Qy_D,U_D,g,nx,ny)
    synchronize()
end