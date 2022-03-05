##############################
#  PLOTTING SCRIPTS
##############################

using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")
np = pyimport("numpy")

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


# Moment calculation
function moment()
    plot_params()

end

# Plot alpha and Vfmax on the same plot
function healing_analysis(Vf, alphaa, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vf, lw = 2.0, label="Max. Slip rate")
    lab1 = "Max. slip rate"
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #ax.set_xlim([0, 250])
    
    col="tab:red"
    ax2 = ax.twinx()
    
    ax2.plot(t./yr2sec, alphaa.*100, c=col, lw=2.0, label="Shear modulus ratio")
    lab2 = "Shear modulus ratio"
    ax.set_xlabel("Time (years)")
    ax2.set_ylabel("Shear Modulus (% of host rock)")
    ax2.set_ylim([99.4, 100.1])
    #  ax2.set_ylim([75, 100])
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)

    #  ax.legend([lab1, lab2], loc=0)
    show()
    
    figname = string(path, "healing_analysis.png")
    fig.savefig(figname, dpi = 300)
end

# Plot the stressdrops after each earthquake
function stressdrop_2(taubefore, tauafter, FltX)
    plot_params()
  
    i = 1;
    #  for i in 1:length(stressdrops[1,:])
      fig = PyPlot.figure(figsize=(7.2, 4.45));
      ax = fig.add_subplot(111);
      ax.plot(taubefore, FltX, lw = 2.0, color="tab:orange", 
              label="Shear stress before the earthquake", alpha=1.0);
      ax.plot(tauafter, FltX, lw = 2.0, color="tab:blue", 
              label="Shear stress after the earthquake", alpha=1.0);
      ax.set_xlabel("Stress drop (MPa)");
      ax.set_ylabel("Depth (km)");
      ax.set_ylim([0,24]);
      ax.set_xlim([15,45]);
      ax.invert_yaxis();
      plt.legend();
      show()
      
      figname = string(path, "shear_stress_im_",i,".png");
      fig.savefig(figname, dpi = 300);
    #  end
end


# Plot shear stress comparison
function shear_stress_comp(shear1b, shear1a, shear2b, shear2a, FltX1, FltX2)
    plot_params()
   
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    ax.plot(shear1b, FltX1, lw = 2.0, color="tab:blue",ls=:dashed, 
            label="Immature Fault Zone: before")
    ax.plot(shear1a, FltX1, lw = 2.0, color="tab:blue", label="Immature Fault Zone: after")
    ax.plot(shear2b, FltX2, lw = 2.0, color="tab:orange", ls=:dashed, 
            label="Mature Fault Zone: before")
    ax.plot(shear2a, FltX2, lw = 2.0, color="tab:orange", label="Mature Fault Zone: after")
    ax.set_xlabel("Shear stress (MPa)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    ax.invert_yaxis()
    plt.legend()
    show()
    
    figname = string(path, "Shear_Stress_002.png");
    fig.savefig(figname, dpi = 300);

end

# Plot rupture_length vs event time
function stem_plot(rl1, rl2, rl3, rl4)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.stem(-rl1./1e3, linefmt="C0-", markerfmt="C0o", basefmt=:None, label="Healing time = 10 yr")
    ax.stem(-rl2./1e3, linefmt="C1-", markerfmt="C1o", basefmt=:None, label="Healing time = 12 yr")
    ax.stem(-rl3./1e3, linefmt="C2-", markerfmt="C2o", basefmt=:None, label="Healing time = 15 yr")
    #  ax.stem(-rl4./1e3, linefmt="C3-", markerfmt="C3o", basefmt=:None, label="Healing time = 20 years")
    ax.set_xlabel("Event number")
    #  ax.set_xlabel("Event time (years)")
    ax.set_ylabel("Rupture length (km)")
    #  ax.set_xlim([10,300])
    plt.legend()
    show()
    
    figname = string(path, "stem_plot2.png")
    fig.savefig(figname, dpi = 300)

end

# Plot slip vs event number
function slipPlot(delfafter2, rupture_len, FltX, Mw, tStart)
    plot_params()
    fig, ax = PyPlot.subplots(nrows=1, ncols=4, sharex="all", sharey="all", figsize=(9.2, 5.00))

    xaxis = tStart[Mw .>2.0]
    delfafter = delfafter2[:,Mw .> 2.0]
    Mw2 = Mw[Mw .> 2.0]

    # Normalize colorbar
    norm = matplotlib.colors.Normalize(vmin = minimum(Mw2), 
                                       vmax=maximum(Mw2)) 
    colors = matplotlib.cm.inferno_r(norm(Mw2))

    ax[1].barh(xaxis, delfafter[end-1,:], height=6, 
              color=colors, align="center"); 
    ax[1].set_ylabel("Time (yr)")
    ax[1].invert_yaxis()
    ax[1].set_title("At 60 m depth")

    trench_depth1 = findall(abs.(FltX) .< 4.0e3)[1]
    trench_depth2 = findall(abs.(FltX) .< 6.0e3)[1]
    trench_depth3 = findall(abs.(FltX) .< 8.0e3)[1]
    
    ax[2].barh(xaxis, delfafter[trench_depth1,:], height=6, 
              color=colors, align="center"); 
    ax[2].invert_yaxis()
    ax[2].set_title("At 4 km depth")
    
    ax[3].barh(xaxis, delfafter[trench_depth2,:], height=6, 
              color=colors, align="center"); 
    ax[3].invert_yaxis()
    ax[3].set_title("At 6 km depth")
    
    ax[4].barh(xaxis, delfafter[trench_depth3,:], height=6, 
              color=colors, align="center"); 
    ax[4].invert_yaxis()
    ax[4].set_title("At 8 km depth")
    
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap="inferno_r")
    sm.set_array([])
    fig.colorbar(sm, shrink=0.9, label="Mw")
    plt.xlabel("Coseismic Slip (m)")
    plt.tight_layout()
    show()
    
    figname = string(path, "coseismic_slip.png")
    fig.savefig(figname, dpi = 300)
end

function plotHypo(hypo)  #S, Slip, SlipVel, Stress, time_)

    # Plot hypocenter
    hist = fit(Histogram, hypo./1e3, closed=:right, nbins=10)

    fig = PyPlot.figure(figsize=(5,7))
    ax = fig.add_subplot(111)

    ax.barh(hist.edges[1][1:end-1], hist.weights)
    #  ax[:plot](collect(1:80), -8*ones(80), "k--", label="Fault Zone Depth")
    ax.set_xlabel("Number of Earthquakes")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Hypocenter Location")
    ax.set_ylim([0, 24])
    ax.invert_yaxis()
    fig.tight_layout()
    show()

    figname = string(path, "hypo.png")
    fig.savefig(figname, dpi = 300)

end

#...........
# Plot MFD
#...........
function MwPlot(Mw)

    hist = fit(Histogram, Mw, nbins = 20)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(8,7))
    ax = fig.add_subplot(111)

    #  ax.plot](hist.edges[1][1:end-1], hist.weights, ".", label="Non-cumulative")
    ax.plot(hist.edges[1][1:end-1], cum, "k.--", markersize=20) #, label="Cumulative")
    ax.set_xlabel("Moment Magnitude (Mw)")
    ax.set_ylabel("Number of Earthquakes")
    ax.set_yscale("log")
    #  ax.set_title("Magnitude-frequency distribution")
    #  ax.set_xlim([2, 7])
    #ax.set_ylim([1, 200])
    #  ax.legend(loc="upper right")
    show()

    figname = string(path, "mfd.png")
    fig.savefig(figname, dpi = 300)
end


# Cumulative sliprate plot
function eqCyclePlot(sliprate, time_stress, FltX)
    indx = findall(abs.(FltX) .<= 19)[1]
    value = sliprate[indx:end,:]
    
    depth = FltX[indx:end]

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    c = ax.imshow(value, cmap="inferno", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e1),
                  interpolation="gaussian",
                  extent=[0,length(sliprate[1,:])/10, 0,19]) 
    
    # for stress
#=       c = ax.imshow(value, cmap="viridis", aspect="auto",
                  vmin=22.5, vmax=40,
                  interpolation="nearest",
                  extent=[0,length(sliprate[1,:])/10, 0,19]) =#
    
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c)
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    show()
    figname = string(path, "mature_stress_2.png")
    fig.savefig(figname, dpi = 300)
    
end


function slip_depth(delfafter, FltX)

    delfafter2 = delfafter[:,481:end]
    FltX2 = FltX[481:end]
    avg_slip = zeros(size(delfafter2)[1], Int(round(481/20)))
    avg_depth = zeros(Int(round(481/20)))
    j = 0
    for i in 1:20:length(FltX2)-1
        j = j+1
        avg_slip[:,j] = mean(delfafter2[:,i:i+20], dims=2)
        avg_depth[j] = mean(FltX2[i:i+20])
    end
    
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)

    #ax.hlines(avg_slip[4,:], 0, avg_depth)
    ax.hlines(avg_depth, 0, avg_slip[4,:], color="tab:blue", lw=2)
    ax.plot(avg_slip[4,:], avg_depth, "D")  # Stem ends
    ax.set_ylim([0,20])
    ax.invert_yaxis()
    ax.set_xlabel("Coseismic Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Earthquake 4")
    show()

    figname = string(path, "slip_depth_4.png")
    fig.savefig(figname, dpi = 300)

end

function x_log_plot(Vf1s, t1s)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)

    t1s_ = abs.(t1s .- t1s[end])
    
    #ax.plot(time_, Vfmax, color="blue", alpha=0.8)
    ax.plot(t1s_, Vf1s, color="tab:blue", label="1 sec", alpha=1.0)
    #ax.plot(t1d_, Vf1d, color="tab:green", label="1 day", alpha=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Max. Slip Rate (m/s)")
    ax.set_xlabel("Time from mainshock (s)")
    #ax.set_xlim([1e9,1e-2])
    
    ax.invert_xaxis()
    x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
    ax.xaxis.set_major_locator(x_major)
    x_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = collect(1:10)*0.1, numticks = 10)
    ax.xaxis.set_minor_locator(x_minor)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    show()
    figname = string(path, "sliprate_log.png")
    fig.savefig(figname, dpi = 300)
    
end

function x_log_plot_2(Vf1s, t1s, Vf5hr, t5hr, Vfno, tno, Vf20d, t20d)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)

    t1s_ = abs.(t1s .- t1s[end])
    t5hr_ = abs.(t5hr .- t5hr[end])
    tno_ = abs.(tno .- tno[end])
    t20d_ = abs.(t20d .- t20d[end])
    
    ax.plot(t1s_, Vf1s, color="tab:blue", label="1 sec", alpha=1.0)
    ax.plot(t5hr_, Vf5hr, color="tab:red", label="5 hr", alpha=1.0)
    ax.plot(tno_, Vfno, color="tab:green", label="No Precursor", alpha=1.0)
    ax.plot(t20d_, Vf20d, color="black", label="10 days", alpha=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Max. Slip Rate (m/s)")
    ax.set_xlabel("Time from mainshock (s)")
    #ax.set_xlim([1e9,1e-2])
    
    ax.invert_xaxis()
    x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
    ax.xaxis.set_major_locator(x_major)
    x_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = collect(1:10)*0.1, numticks = 10)
    ax.xaxis.set_minor_locator(x_minor)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    plt.legend()
    show()
    figname = string(path, "sliprate_log.png")
    fig.savefig(figname, dpi = 300)
    
end

function contourf_plot(sliprate, time_stress, FltX)
    indx = findall(abs.(FltX) .<= 20)[1]
    value = sliprate[:, indx:end]

    time2 = abs.(time_stress .- time_stress[end])
    
    depth = FltX[indx:end]

    X, Y = np.meshgrid(time2, depth) 

    #return X, Y, value
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    c = ax.contourf(X, Y, abs.(value)', cmap="inferno", levels=10 .^(LinRange(-12, 2, 50)),
                    norm=matplotlib.colors.LogNorm(), extend="max")
    
    ax.set_xlabel("Time from mainshock (s)")
    
    ax.set_ylabel("Depth (km)")
    #ax.set_ylim([5,15])

    #ax.set_xscale("log")
    #ax.invert_xaxis()
#=     x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
    ax.xaxis.set_major_locator(x_major)
    x_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = collect(1:10)*0.1, numticks = 10)
    ax.xaxis.set_minor_locator(x_minor)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()) =#

    ax.invert_yaxis()
    # cbaxes = fig.add_axes([5, 0.15, 0.7, 0.03])
    plt.colorbar(c, label="Slip Rate (m/s)", ax=ax)
    #cbar = fig.colorbar(c)
    # plt.tight_layout()
    #cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    show()
    figname = string(path, "sliprate_evno_004.png")
    fig.savefig(figname, dpi = 300)
    
end


# Plot Vfmax
function VfmaxPlot_solo(Vfmax, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vfmax, color="tab:blue", label="No Precursor")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_xlim([0,500])
    show()
    
    figname = string(path, "Vfmax03.png")
    fig.savefig(figname, dpi = 300)
end
function VfmaxPlot(Vfmax, t, Vf2, t2, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vfmax, "r", lw = 2.0, label="No Precursor")
    ax.plot(t2./yr2sec, Vf2, "b--", lw = 2.0, label="Precursor 1 sec. before EQ.")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_xlim([0,500])
    plt.legend(loc="upper right")
    show()
    
    figname = string(path, "Vfmax03.png")
    fig.savefig(figname, dpi = 300)
end

# Plot Vsurface
function VsurfPlot(Vsurf10, Vsurf12, Vsurf15, t10, t12, t15, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(8.2, 6.00))
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    ax.plot(t10./yr2sec, Vsurf10, lw = 2.0, label="Healing time = 10 yr")
    ax2.plot(t12./yr2sec, Vsurf12, lw = 2.0, label="Healing time = 12 yr")
    ax3.plot(t15./yr2sec, Vsurf15, lw = 2.0, label="Healing time = 15 yr")
    ax3.set_xlabel("Time (years)")
    ax2.set_ylabel("Surface. Slip rate (m/s)")
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax3.set_yscale("log")
    #  ax.set_xlim([230,400])
    #  plt.legend()
    plt.tight_layout()
    show()
    
    figname = string(path, "Vsurf02.png")
    fig.savefig(figname, dpi = 300)
end

# Plot seismic moment
function moment_plot(Mo_a1, Mo_a2, Mo_a3, Mo_s1, Mo_s2, Mo_s3)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(Mo_a1./1e0, "bo", label="Aseismic Moment: no precursor")
    ax.plot(Mo_a2./1e0, "ro", label="Aseismic Moment: 0.1% precursor")
    ax.plot(Mo_a3./1e0, "go", label="Aseismic Moment: 0.5% precursor")
    #ax.plot(Mo_s1./1e0, "bo", alpha=0.4, label="Seismic Moment: no precursor")
    #ax.plot(Mo_s2./1e0, "ro", alpha=0.4,label="Seismic Moment: 0.1% precursor")
    #ax.plot(Mo_s3./1e0, "go", alpha=0.4,label="Seismic Moment: 0.5% precursor")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Moment")

    plt.legend()
    show()
    figname = string(path, "amoment_02.png")
    fig.savefig(figname, dpi = 300)

end

# Plot alpha
function alphaaPlot(alphaa, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(t./yr2sec, alphaa.*100, lw = 2)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Contrast (%)")
    #  ax.set_xlim([230,400])
    show()


    figname = string(path, "alpha_01.png")
    fig.savefig(figname, dpi = 300)
end

# Compare alpha
function VfComp(Vfmax1, Vfmax2, Vfmax3, t1, t2, t3, tS1, tS2, tS3, tE1, tE2, tE3, yr2sec)
    plot_params()


	idS1 = findmax(t1 .>= tS1[6])[2]
	idE1 = findmax(t1 .>= tE1[6])[2]
	
	# 1 sec precursor
	idS2 = findmax(t2 .>= tS2[9])[2]
	idE2 = findmax(t2 .>= tE2[9])[2]	

    # 12 days precursor
    idS3 = findmax(t3 .>= tS3[9])[2]
	idE3 = findmax(t3 .>= tE3[9])[2]	

	# ref is the translation factor to bring both plots to
	# the same zero (in seconds)
	ref1 = -30
	#ref1 = -0
	ref2 = -23.5

    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t1[idS1-300:idE1] .- t1[idS1], Vfmax1[idS1-300:idE1], lw = 2.0, label="No Precursor")
    ax.set_ylabel("Max. Slip rate (m/s)")
    # ax.set_xlim([0, 250])
    
    #col="tab:red"
    #ax2 = ax.twinx()
    
    ax.plot(t2[idS2-450:idE2] .- t2[idS2] .- ref1, Vfmax2[idS2-450:idE2], lw=2.0, alpha=0.8, label="Precursor onset = 1 sec")
    ax.plot(t3[idS3-150:idE3] .- t3[idS3] .- ref1, Vfmax3[idS3-150:idE3], lw=2.0, alpha=0.8, label="Precursor onset = 12 days")
    ax.set_xlabel("Time (sec)")
    ax.set_yscale("log")
    # ax2.set_ylabel("Max. acceleration (m/s2)")
    #ax2.set_yscale("log")
    # ax2.set_ylim([35, 60])
    #  ax2.set_ylim([75, 100])
    #ax2.get_xaxis().set_tick_params(color=col)
    #ax2.tick_params(axis="x", labelcolor=col)

    ax.legend(loc=0)
    show()

    figname = string(path, "aseismic_Vf_004.png")
    fig.savefig(figname, dpi = 300)
end

# Plot cumulative slip
function cumSlipPlot(delfsec, delfyr, FltX)
    indx = findall(abs.(FltX) .<= 18)[1]

    delfsec2 = transpose(delfsec[:,indx:end])
    delfyr2 = transpose(delfyr)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    plt.rc("font",size=12)

    ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
    ax.plot(delfsec2, FltX[indx:end], color="chocolate", lw=1.0)
    ax.set_xlabel("Accumulated Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    ax.set_xlim([0,maximum(delfsec2)])
    ax.set_xlim([0,12.0])
    
    ax.invert_yaxis()
    
    show()
    
    figname = string(path, "cumulative_slip.png")
    fig.savefig(figname, dpi = 300)

end

# Plot friction parameters
function icsPlot(a_b, Seff, tauo, FltX)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(Seff, FltX, "k-", label="Normal Stress")
    ax.plot(tauo, FltX, "k--", label="Shear Stress")
    ax.set_xlabel("Stresses (MPa)")
    ax.set_ylabel("Depth (km)")
    plt.legend(loc="lower right") 
    
    col="tab:blue"
    ax2 = ax.twiny()
    ax2.plot(a_b, FltX, label="(a-b)")
    ax2.set_xlabel("Rate-state friction value", color=col)
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)
    
    ax.set_ylim([0,16])
    ax.invert_yaxis()
    show()
    
    figname = string(path, "ics_02.png")
    fig.savefig(figname, dpi = 300)
end


# Plot stressdrop comparison
function sd_comp(sd8, sd9, sd10, sd11, FltX)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(sd8, FltX, color="k", lw=0.4, label="1 sec before EQ")
    ax.plot(sd9, FltX, color="blue", lw=0.4, label="5 hrs before EQ")
    ax.plot(sd10, FltX, color="red", lw=0.4, label="1 day before EQ")
    ax.plot(sd11, FltX, color="g", lw=0.2, label="20 days before EQ")
    ax.set_xlabel("Shear Stress (MPa)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,16])
    ax.invert_yaxis()
    plt.legend() 
    show()
    
    figname = string(path, "shear_stress_eq_02.png")
    fig.savefig(figname, dpi = 300)
end
