# initialisation & global definition(s)
using Plots, LaTeXStrings, ProgressMeter, ColorSchemes, DelimitedFiles, CSV

# include dependencies & function call(s)
include("../src/fun/hillshade.jl")
include("../src/fun/plots.jl")

@views function main()
    # physical constant
    z = readdlm("./docs/example/hillshade/data/dtm_1m.txt")
    nx,ny = size(z)
    Δx = 1.0
    Δy = Δx
    xc = 1:Δx:nx*Δx
    yc = 1:Δx:ny*Δy
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
    @time hs=hillshade(xc,yc,z,Δx,Δy,45.0,315.0,nx,ny)
    @time hillshade_plot(yc,xc,hs',45.0,315.0,0.75)
    savefig("./src/out/hillshade_from_DEM.png")



    
end
main()
# https://techytok.com/lesson-parallel-computing
# https://nbviewer.org/github/daniel-koehn/Differential-equations-earth-system/blob/master/10_Shallow_Water_Equation_2D/01_2D_Shallow_Water_Equations.ipynb