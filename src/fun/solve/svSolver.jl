# global variable(s) & global declaration(s)
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=false
    )
ϵ    = 1.0e-10
path = "src/out/"

# include dependencies & function call(s) for svSolver.jl
include("../plots.jl")
include("../hillshade.jl")
include("advSolve.jl")
include("souSolve.jl")
include("get.jl")

@views function svSolver(xc,yc,h,Qx,Qy,z,g,CFL,T,tC,Δx,Δy,nx,ny,Dsim)
    solv_type  = Dsim.solv_type
    make_gif   = Dsim.make_gif
    flow_type  = Dsim.flow_type
    pcpt_onoff = Dsim.pcpt_onoff
    println("[=> generating initial plots & exporting...")
    # display initial stuffs
    η0   = minimum(h.+z)
    zmin = minimum(z)
    ηmax0= maximum(h.+z)
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        z_plot(xc,yc,z)
    savefig(path*"plot_z_init.png")
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        h_plot(xc,yc,h,maximum(h),nx,ny,0.0,flow_type)
    savefig(path*"plot_h_init.png")
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        free_surface_plot(xc,yc,h,z,η0,0.75*(ηmax0-η0),nx,ny,0.0)
    savefig(path*"plot_eta_init.png")
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        wave_plot(xc,yc,h,z,η0,(ηmax0-η0),nx,ny,0.0)
    savefig(path*"plot_wave_height_init.png")
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        profile_plot(xc,yc,h,z,zmin,10.0,nx,ny,0.0)
    savefig(path*"plot_profile_init.png")   
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
        hs=hillshade(xc,yc,z,Δx,Δy,45.0,315.0,nx,ny)
        hillshade_plot(xc,yc,hs,45.0,315.0,0.75)
    savefig(path*"plot_hillshade.png")

    # set & get vectors
    U,F,G = getUF(h,Qx,Qy,g,nx,ny)
    # set time
    t     = 0.0
    # plot & time stepping parameters
    it    = 0
    ctr   = 0
    # generate GIF
    if make_gif==true
        println("[=> initializing & configuring .gif...")
        anim = Animation()
    end
    # action
    println("[=> action!")
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while t<T
    	# adaptative Δt
        Δt  = getΔt(h,Qx,Qy,g,Δx,Δy,CFL,nx,ny)
        # advection step solution
        h,Qx,Qy = advSolve(h,Qx,Qy,z,U,F,G,g,Δx,Δy,Δt,nx,ny,solv_type)
        # source step solution
        h,Qx,Qy = souSolve(h,Qx,Qy,z,U,F,G,g,Δx,Δy,t,Δt,nx,ny,flow_type,pcpt_onoff)
        # update current time
        t  += Δt
        it += 1
        if t > ctr*tC
            ctr+=1
                fig=gr(size=(2*250,2*125),markersize=2.5)       
                    #fig=wave_plot(xc,yc,h,z,η0,0.1*ηmax0,nx,ny,t)
                    #fig=free_surface_plot(xc,yc,h,z,η0,0.25*(maximum(h.+z)-η0),nx,ny,t)
                    #fig=discharge_plot(xc,yc,h,Qx,Qy,z,2.5,nx,ny,t)
                    #fig=profile_plot(xc,yc,h,z,zmin,10.0,nx,ny,t)
                    fig=h_plot(xc,yc,h,0.5,nx,ny,t,flow_type)
                    if make_gif==true
                        frame(anim,fig)
                    end
        end
        next!(prog;showvalues = [("[nx,ny]",(nx,ny)),("iteration(s)",it),("(✗) t/T",round(t/T,digits=2))])
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = [("[nx,ny]",(nx,ny)),("iteration(s)",it),("(✓) t/T",1.0)])
    println("[=> generating final plots, exporting & exiting...")
    if make_gif==true
        gif(anim,path*solv_type*"_wave.gif")
    end
    savefig(path*solv_type*"_plot.png")

    free_surface_plot(xc,yc,h,z,η0,0.3*(maximum(h.+z)-η0),nx,ny,t)
    savefig(path*solv_type*"_freesurface.png")

end