using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")
using3D()

# from mpl_toolkits.mplot3d import Axes3D

# Default plot params
function plot_params()
  plt.rc("xtick", labelsize=16)
  plt.rc("ytick", labelsize=16)
  plt.rc("xtick", direction="in")
  plt.rc("ytick", direction="in")
  plt.rc("font", size=15)
  plt.rc("figure", autolayout="True")
  plt.rc("axes", titlesize=16)
  plt.rc("axes", labelsize=17)
  plt.rc("xtick.major", width=1.5)
  plt.rc("xtick.major", size=5)
  plt.rc("ytick.major", width=1.5)
  plt.rc("ytick.major", size=5)
  plt.rc("lines", linewidth=2)
  plt.rc("axes", linewidth=1.5)
  plt.rc("legend", fontsize=13)
  plt.rc("mathtext", fontset="stix")
  plt.rc("font", family="STIXGeneral")

  # Default width for Nature is 7.2 inches, 
  # height can be anything
  #plt.rc("figure", figsize=(7.2, 4.5))
end

# Plot shear stress comparison
function Vfmax_comp(Vfmax1, Vfmax2, t1, t2, tS1, tS2, tE1, tE2)
    plot_params()
 
	# Without vs change
	idS1 = findmax(t1 .>= tS1[5])[2]
	idE1 = findmax(t1 .>= tE1[5])[2]
	
	# With vs change
	idS2 = findmax(t2 .>= tS2[2])[2]
	idE2 = findmax(t2 .>= tE2[2])[2]	

	# ref is the translation factor to bring both plots to
	# the same zero (in seconds)
	ref = -23.5

    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    ax.plot(t1[idS1:idE1] .- t1[idS1] .- ref, Vfmax1[idS1:idE1], lw = 2.0, color="tab:blue", 
            label="1 km DFZ, 0% vs contrast")
    ax.plot(t2[idS2:idE2] .- t2[idS2], Vfmax2[idS2:idE2], lw = 2.0, color="tab:orange", alpha = 0.9,
            label="1 km DFZ, 1% vs contrast") 
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Max. slip rate (m/s)")
    plt.legend()
    show()
    
    figname = string(path, "VfComp_01.png");
    fig.savefig(figname, dpi = 300);

end


# Plot shear stress comparison
function Vfmax_comp_full(t1, Vf1, t2, Vf2, t3, Vf3, t4, Vf4, t5, Vf5, yr2sec)
    plot_params()
 
    fig = PyPlot.figure(figsize=(7.2, 3.00))
    ax = fig.add_subplot(111)
#=     ax.plot(t1, Vf1, lw = 2.0, color="black", 
            label="1 sec. before EQ") =#
    ax.plot(t2, Vf2, lw = 2.0, color="blue", 
            label="5 hrs before EQ")
   ax.plot(t3, Vf3, lw = 2.0, color="red",
            label="1 day before EQ") 
   ax.plot(t4, Vf4, lw = 2.0, color="g",
            label="20 days before EQ") 
#=    ax.plot(t5, Vf5, lw = 2.0, color="indigo",
            label="45 days before EQ")  =#
    ax.set_xlabel("Time (yrs)")
    ax.set_ylabel("Max. slip rate (m/s)")
	ax.set_yscale("log")
    ax.set_xlim([200, 225])
    # ax.set_xlim([0, 400])
    # plt.legend()
    show()
    
    figname = string(path, "Vfmax_comp_full_6.png");
    fig.savefig(figname, dpi = 300);

end

# Try to plot rupture speed
function rupture_speed(delfsec, FltX)
    plot_params()

    # Rupture speed calculation
    indx = findall(abs.(FltX) .<= 18)[1]
    delfsec2 = transpose(delfsec[:,indx:end])
    FltX2 = FltX[indx:end]
    len = size(delfsec2)[2]

    slip = zeros(len,2)
    depth = zeros(len,2)
    id1 = 0; id2 = 0
    for i in 1:len
        id1 = findall(delfsec2[:,i] .> 0.95*maximum(delfsec2[:,i]))[1]
        id2 = findall(delfsec2[:,i] .> 0.95*maximum(delfsec2[:,i]))[end]

        slip[i,1] = delfsec2[id1,i]
        slip[i,2] = delfsec2[id2,i]
        
        depth[i,1] = FltX2[id1]
        depth[i,2] = FltX2[id2]
    end
    
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)

    ax.scatter(slip[:,1], depth[:,1], s=4, color="grey")
    ax.scatter(slip[:,2], depth[:,2], s=4, color="grey")
    #ax.plot(slip[:,1], depth[:,1], "-", lw=0.4, color="tab:blue")
    #ax.plot(slip[:,2], depth[:,2], "-", lw=0.4, color="tab:blue")
    ax.plot(delfsec2, FltX2, color="tab:orange", alpha=1.0, lw=0.2)

    ax.set_xlabel("Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_xlim([0,12])
    ax.invert_yaxis()

    figname = string(path, "rupture_speed_2.png");
    fig.savefig(figname, dpi = 300);
end

function surface(time_stress, FltX, sliprate, yr2sec)
    plot_params()
 
    X, Y = np.meshgrid(time_stress, FltX)
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    ax = fig.gca(projection="3d")

    ax.surf(X, Y, sliprate)
    #ax.set_xlabel("Slip (m)")
    #ax.set_ylabel("Depth (km)")
    #ax.set_xlim([0,12])
    #ax.invert_yaxis()

    figname = string(path, "surface.png");
    fig.savefig(figname, dpi = 300);
end

# Plot cumulative slip
function cumSlipPlot2(delfsec, delfyr, FltX)
    indx = findall(abs.(FltX) .<= 18)[1]

    delfsec2 = transpose(delfsec[:,indx:end])
    delfyr2 = transpose(delfyr)

    plot_params()
    fig = PyPlot.figure(figsize=(2.2, 4.45))
    ax = fig.add_subplot(111)
    plt.rc("font",size=12)

    ax.plot(delfyr2, FltX, color="royalblue", lw=0.4)
    ax.plot(delfsec2, FltX[indx:end], color="chocolate", lw=0.4)
    #ax.set_xlabel("Accumulated Slip (m)")
    #ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    ax.set_xlim([0,maximum(delfsec2)])
    ax.set_xlim([7,12.0])
    
    ax.invert_yaxis()
    
    show()
    
    figname = string(path, "temp_cum_slip.png")
    fig.savefig(figname, dpi = 300)

end

# Plot sliprate nucleation
function slipr_nucl(sliprate, FltX)
    plot_params()
    
    indx = findall(abs.(FltX) .<= 24)[1]
    value = transpose(abs.(sliprate[:,indx:end]))

    fig = PyPlot.figure(figsize=(3.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(value, FltX[indx:end], color="k", lw=0.4)
    ax.set_xlabel("Sliprate (m/s) ")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    ax.set_xscale("log")
    ax.invert_yaxis()
    # plt.legend() 
    show()
    
    figname = string(path, "slipr_nucl.png")
    fig.savefig(figname, dpi = 300)
end