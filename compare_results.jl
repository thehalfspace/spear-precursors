using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")
include("$(@__DIR__)/post/compare_plots.jl")
include("$(@__DIR__)/post/vs_change_scripts.jl")
include("$(@__DIR__)/post/output_seismograms.jl")

# path to save files
global path = "$(@__DIR__)/plots/vs_damage_comp/"
mkpath(path)

global out_path4 = "$(@__DIR__)/data/vs_damage_04/"
global out_path8 = "$(@__DIR__)/data/vs_damage_08/"
global out_path9 = "$(@__DIR__)/data/vs_damage_09/"
global out_path10 = "$(@__DIR__)/data/vs_damage_10/"
global out_path11 = "$(@__DIR__)/data/vs_damage_11/"

# Global variables
yr2sec = 365*24*60*60

# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path4, "params.out"), header=false)

Seff = params[1,:]
tauo = params[2,:]
FltX = params[3,:]
cca = params[4,:]
ccb = params[5,:]
Lc = params[6,:]

# Index of fault from 0 to 18 km
flt18k = findall(FltX .<= 18)[1]

time_vel4 = readdlm(string(out_path4, "time_velocity.out"), header=false)
t4 = time_vel4[:,1]
Vfmax4 = time_vel4[:,2]
Vsurface4 = time_vel4[:,3]
alphaa4 = time_vel4[:,4]

time_vel8 = readdlm(string(out_path8, "time_velocity.out"), header=false)
t8 = time_vel8[:,1]
Vfmax8 = time_vel8[:,2]
Vsurface8 = time_vel8[:,3]
alphaa8 = time_vel8[:,4]

time_vel9 = readdlm(string(out_path9, "time_velocity.out"), header=false)
t9 = time_vel9[:,1]
Vfmax9 = time_vel9[:,2]
Vsurface9 = time_vel9[:,3]
alphaa9 = time_vel9[:,4]

time_vel10 = readdlm(string(out_path10, "time_velocity.out"), header=false)
t10 = time_vel10[:,1]
Vfmax10 = time_vel10[:,2]
Vsurface10 = time_vel10[:,3]
alphaa10 = time_vel10[:,4]

time_vel11 = readdlm(string(out_path11, "time_velocity.out"), header=false)
t11 = time_vel11[:,1]
Vfmax11 = time_vel11[:,2]
Vsurface11 = time_vel11[:,3]
alphaa11 = time_vel11[:,4] 

event_time4 = readdlm(string(out_path4, "event_time.out"), header=false)
tStart4 = event_time4[:,1]
tEnd4 = event_time4[:,2]
hypo4 = event_time4[:,3]

event_time8 = readdlm(string(out_path8, "event_time.out"), header=false)
tStart8 = event_time8[:,1]
tEnd8 = event_time8[:,2]
hypo8 = event_time8[:,3]

event_time9 = readdlm(string(out_path9, "event_time.out"), header=false)
tStart9 = event_time9[:,1]
tEnd9 = event_time9[:,2]
hypo9 = event_time9[:,3]

event_time10 = readdlm(string(out_path10, "event_time.out"), header=false)
tStart10 = event_time10[:,1]
tEnd10 = event_time10[:,2]
hypo10 = event_time10[:,3]

event_time11 = readdlm(string(out_path11, "event_time.out"), header=false)
tStart11 = event_time11[:,1]
tEnd11 = event_time11[:,2]
hypo11 = event_time11[:,3]