# include dependencies & function call(s)
include("../bc/getBCs.jl")
include("../flux/fluxes.jl")
include("../upd/update.jl")
@views function advSolve(h,Qx,Qy,z,U,F,G,g,Δx,Δy,Δt,nx,ny,type)
    # assembly of conservative variables vector and flux function vector
        getU!(U,h,Qx,Qy,nx,ny)
    # ghost cells
        #2) periodic in x
        Ubc = getBCs(U,nx,ny,3,"outflow")
        zbc = getBCs(z,nx,ny,1,"dirichlet")          
    # find inter-cell fluxes (e.g., HLL-type fluxes or the Rusanov fluxes)
        fluxes!(F,Ubc,zbc,g,Δx,Δt,type,nx,ny,"x")
    # update conservative variable vector, e.g., splitting scheme
        updateU!(U,F,(Δt/Δx),nx,ny,3)
    # ghost cells    
        #2) periodic in y
        Ubc = getBCs(permutedims(U,(2,1,3)),ny,nx,3,"outflow")
        zbc = getBCs(permutedims(z,(2,1)  ),ny,nx,1,"dirichlet")
    # find inter-cell fluxes (e.g., HLL-type fluxes or the Rusanov fluxes)
        fluxes!(G,Ubc,zbc,g,Δy,Δt,type,ny,nx,"y")
    # update conservative variable vector, e.g., splitting scheme
        updateU!(U,permutedims(G,(2,1,3)),(Δt/Δy),nx,ny,3)
    return copy(U[:,:,1]),copy(U[:,:,2]),copy(U[:,:,3])
end