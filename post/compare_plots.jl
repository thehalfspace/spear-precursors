##############################
#  PLOTTING SCRIPTS
##############################

using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")

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

# Plot Vfmax
function VfmaxComp(Vf1, t1, Vf2, t2, Vf3, t3, Vf4, t4, Vf5, t5, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t1./yr2sec, Vf1, lw = 2.0, label="Homogeneous, 0% red.")
    ax.plot(t2./yr2sec, Vf2, lw = 2.0, label="Homogeneous, 0.1% red.")
    # ax.plot(t3./yr2sec, Vf3, lw = 2.0, label="inf/0% DFZ, 1% red.")
    # ax.plot(t4./yr2sec, Vf4, lw = 2.0, label="1 km DFZ, 1% red.")
    # ax.plot(t5./yr2sec, Vf5, lw = 2.0, label="1 km DFZ, 2% red.")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #  ax.set_xlim([230,400])
	ax.legend(loc="upper left")
    show()
    
    figname = string(path, "Vfcomp09.png")
    fig.savefig(figname, dpi = 300)
end

function VfEvent(Vf1, t1, Vf2, t2, Vf3, t3, Vf4, t4, Vf5, t5, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t1./yr2sec, Vf1, lw = 2.0, label="Homogeneous")
    ax.plot(t2./yr2sec, Vf2, lw = 2.0, label="2 km DFZ, 1% red.")
    # ax.plot(t3./yr2sec, Vf3, lw = 2.0, label="2 km DFZ, 2% red.")
    # ax.plot(t4./yr2sec, Vf4, lw = 2.0, label="1 km DFZ, 1% red.")
    # ax.plot(t5./yr2sec, Vf5, lw = 2.0, label="1 km DFZ, 2% red.")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #  ax.set_xlim([230,400])
	ax.legend()
    show()
    
    figname = string(path, "Vfevent03.png")
    fig.savefig(figname, dpi = 300)
end

function Vfmax_comp(Vfmax1, Vfmax2, t1, t2, tS1, tS2, tE1, tE2)
    plot_params()
 
	# Without vs change
	idS1 = findmax(t1 .>= tS1[2])[2]
	idE1 = findmax(t1 .>= tE1[2])[2]
	
	# With vs change
	idS2 = findmax(t2 .>= tS2[2])[2]
	idE2 = findmax(t2 .>= tE2[2])[2]	

	# ref is the translation factor to bring both plots to
	# the same zero (in seconds)
	ref = -0.8

    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    ax.plot(t2[idS2:idE2] .- t2[idS2], Vfmax2[idS2:idE2], lw = 2.0, color="tab:orange", alpha = 1.0,
            label="Homogeneous, 0.1% red.") 
    ax.plot(t1[idS1:idE1] .- t1[idS1] .- ref, Vfmax1[idS1:idE1], lw = 2.0, color="tab:blue", alpha = 0.8, 
            label="Homogeneous, 0% red.")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Max. slip rate (m/s)")
    # ax.set_yscale("log")
    plt.legend()
    show()
    
    figname = string(path, "VfComp_05.png");
    fig.savefig(figname, dpi = 300);

end

