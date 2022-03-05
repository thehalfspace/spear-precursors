##############################
#  PLOT OUTPUT SEISMOGRAMS
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

function get_event_index(vel_seis)
    l = length(vel_seis[:,1])
    event_index = zeros(l)

    evno = 1    # event number
    # separate each event
    for i in 1:l
        if vel_seis[i,1] == "end"
            event_index[evno] = i-1
            evno = evno + 1
        end
    end
    return Int.(filter(!iszero, event_index))
end


function vel_seis_plot(vel_seis)
    event_index = get_event_index(vel_seis)

    # first evet
    vel_seis_2nd = vel_seis[1:event_index[1],:]
    #  vel_seis_2nd = vel_seis[event_index[1]+2:event_index[2],:]
   
    # X axis time in seconds
    x = collect(1:length(vel_seis_2nd[:,1]))./10
    
    plot_params()
    fig, ax = plt.subplots(nrows=5, ncols=1, sharex="all", 
                           sharey="all", figsize=(6, 9))
    
    ax[1].plot(x, vel_seis_2nd[:,1], lw = 2.0, label="50 m from fault")
    ax[2].plot(x, vel_seis_2nd[:,2], lw = 2.0, label="100 m from fault")
    ax[3].plot(x, vel_seis_2nd[:,3], lw = 2.0, label="200 m from fault")
    ax[4].plot(x, vel_seis_2nd[:,4], lw = 2.0, label="500 m from fault")
    ax[5].plot(x, vel_seis_2nd[:,5], lw = 2.0, label="1000 m from fault")

    plt.xlabel("Time (s)")
    ax[3].set_ylabel("Velocity (m/s)") 

    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()
    ax[5].legend()
    #  plt.tight_layout()
    show()
    
    figname = string(path, "vel_seis_01.png")
    fig.savefig(figname, dpi = 300)
end
