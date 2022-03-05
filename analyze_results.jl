using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")
include("$(@__DIR__)/post/output_seismograms.jl")
include("$(@__DIR__)/post/vs_change_scripts.jl")

# path to save files
global path = "$(@__DIR__)/plots/paper_3_ref_model/"
mkpath(path)

global out_path = "$(@__DIR__)/data/paper_3_ref_model/"

# Global variables
yr2sec = 365*24*60*60

# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path, "params.out"), header=false)

Seff = params[1,:]
tauo = params[2,:]
FltX = params[3,:]
cca = params[4,:]
ccb = params[5,:]
Lc = params[6,:]

# Read data
event_time = readdlm(string(out_path, "event_time.out"), header=false)
tStart = event_time[:,1]
tEnd = event_time[:,2]
hypo = event_time[:,3]

event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
indx = Int(length(event_stress[1,:])/2)
taubefore = event_stress[:,1:indx]
tauafter = event_stress[:,indx+1:end]

delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)
#  slip = readdlm(string(out_path, "slip.out"), header=false)
sliprate = readdlm(string(out_path, "sliprate.out"), header=false)

# Index of fault from 0 to 18 km
flt18k = findall(FltX .<= 18)[1]

time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
t = time_vel[:,1]
Vfmax = time_vel[:,2]
Vsurface = time_vel[:,3]
alphaa = time_vel[:,4]
amax = time_vel[:,5]
isolver = time_vel[:,6]


rho1 = 2670
vs1 = 3464
rho2 = 2500
vs2 = 0.6*vs1
mu = rho2*vs2^2

delfsec = readdlm(string(out_path, "delfsec.out"))
delfyr = readdlm(string(out_path, "delfyr.out"))
stress = readdlm(string(out_path, "stress.out"), header=false)
time_stress = readdlm(string(out_path, "time_stress.out"), header=false)
slr_max_01 = maximum(sliprate, dims=2)

start_index = get_index(stress', taubefore')
stressdrops = taubefore .- tauafter


Mw, del_sigma, fault_slip, rupture_len, seismic_moment =
        moment_magnitude_new(mu, FltX, delfafter', stressdrops', t);

# aseismic_moment = get_moment(sliprate, time_stress, Mw, rupture_len, mu, FltX)
# OUTPUT SEISMOGRAMS
vel_seis = readdlm(string(out_path, "velocity_seismogram.out"), header=false);
