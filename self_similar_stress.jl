using LinearAlgebra
using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
using FFTW
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
end


function stress_amplitude(k, ac, H)
    # Correlation length: ac = 2π/kc
    kc = 2π/ac
    τ_k = zeros(size(k))
    for i in 1:size(k)[1]
        for j in 1:size(k)[2]
            if k[i,j] <= kc
                τ_k[i,j] = kc^(-(1+H))
            else
                τ_k[i,j] = k[i,j]^(-(1+H))  
            end
        end
    end
    return τ_k
end

function plot_stress_amp()
    plot_params()

    # input values
    x = LinRange(-1,2,20)
    k = 10 .^(x)
    ac = 5e0    # in meters
    H = 1       # Test case


    tau_k = stress_amplitude(k, ac, H)

    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)


    # ax.plot(k, ifft(tau_k))
    ax.plot(k, tau_k)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_yscale("log")
    ax.set_xscale("log")

    show()
end

function norm_alize(x)
    return (x .- mean(x))/maximum(x)
end

function plot_2D_self_similar_stress()
    plot_params()

    N = 500
    # input values
    x = LinRange(0.1,100,N)
    z = LinRange(0.1,100,N)

    kx, kz = np.meshgrid(x, z)


    # random phase values
    ϕ = 100 * rand(length(x), length(x))

    # k = 10 .^(kx.^2 .+ kz.^2 .+ ϕ.^2)
    k = (kx.^2 .+ kz.^2 .+ ϕ.^2) .^ 0.5
    ac = 5e0    # in meters
    H = 1       # Test case

    tau_k = stress_amplitude(k, ac, H)
    tau_x_z = ifft(tau_k)

    r_tau = abs.(tau_x_z)
    # r_tau = reshape(abs.(tau_x_z), length(x), length(x))

    r_tau2 = r_tau./maximum(r_tau)

    # return r_tau

    # return x_real, x_im
    fig = PyPlot.figure(figsize=(9.2, 6.45))
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)


    extent = [0, 1, 0, 1]

    #plt.clf()
    c = ax.imshow(r_tau2, interpolation="none", aspect="auto", 
        norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=1e-1), 
        extent=extent, cmap="inferno")

    ax.set_xlabel("Strike(normalized)")
    ax.set_ylabel("Dip (normalized)")
    cbar=fig.colorbar(c)
    # ax.set_yscale("log")
    # ax.set_xscale("log")

    ax2.plot(r_tau2[5:1:end,Int(length(x)/2)])
    ax2.set_xlabel("Number of points")
    ax2.set_ylabel("Stress (normalized)")

    return r_tau2[5:1:end, Int(length(x)/2)]

    show()
end