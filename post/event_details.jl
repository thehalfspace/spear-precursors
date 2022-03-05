#################################
# MODULE FOR SOME CALCULATIONS
# FROM SIMULATION OUTPUT
#################################
using LinearAlgebra

# get index of start of rupture
function get_index(seismic_stress, taubefore)

    len = length(taubefore[:,1])
    index_start = zeros(Int, length(taubefore[1,:]))
    for i in 1:length(taubefore[1,:])
        temp = zeros(length(seismic_stress[1,:]))
        #  index_start[i] = findall(seismic_stress[len,:] .== taubefore[len,i])[1]
        for j in 1:length(seismic_stress[1,:])
            temp[j] = norm(seismic_stress[:,j] .- taubefore[:,i])
        end

        index_start[i] = findmin(temp)[2]
    end

    index_start
end

function get_moment(sliprate, time_, Mw, rupture_len, mu, FltX)
    slr_max = maximum(sliprate, dims=2)

    aseismic_id = Tuple.(findall( 1e-3 .> slr_max .> 1e-7))
    seismic_id = Tuple.(findall( slr_max .> 1e-3))

    aseismic_inds = [aseismic_id[i][1] for i in 1:length(aseismic_id)]
    seismic_inds = [seismic_id[i][1] for i in 1:length(seismic_id)]

    aseismic_dt = diff(time_[aseismic_inds])
    seismic_dt = diff(time_[seismic_inds])

    aseismic_slr_temp = sliprate[aseismic_inds,:]
    seismic_slr_temp = sliprate[seismic_inds,:]

    aseismic_slr = aseismic_slr_temp[2:end,:]
    seismic_slr = seismic_slr_temp[2:end,:]
    
    aseismic_slip = aseismic_dt.*aseismic_slr
    seismic_slip = seismic_dt.*seismic_slr

    Mo_a_temp = -aseismic_slip[:,2:end]' .* diff(FltX) .* mu

    Mo_a_2 = sum(Mo_a_temp, dims=1)

    # return Mo_a_2

    Mo_a = zeros(4)
    Mo_a[1] = sum(Mo_a_2[1:20])
    Mo_a[2] = sum(Mo_a_2[21:50])
    Mo_a[3] = sum(Mo_a_2[51:150])
    Mo_a[4] = sum(Mo_a_2[151:end])
    
    Mw_a = zeros(4)
#=     Mw_a[1] = (2/3)*log10.(Mo_a1.*1e7) .- 10.7
    Mw_a[2] = (2/3)*log10.(Mo_a2.*1e7) .- 10.7
    Mw_a[3] = (2/3)*log10.(Mo_a3.*1e7) .- 10.7
    Mw_a[4] = (2/3)*log10.(Mo_a4.*1e7) .- 10.7 =#
    # Calculate aseismic moment
    rup_len = rupture_len[Mw .> 4.5]

    return Mo_a
end

#.................................................
# Compute the final Coseismic slip for each event
#.................................................
function Coslip(S, Slip, SlipVel, Stress, time_=zeros(1000000))
    Vfmax = maximum(SlipVel, dims = 1)[:]

    delfafter::Array{Float64,2} = zeros(size(Slip))
    tStart::Array{Float64} = zeros(size(Slip[1,:]))
    tEnd::Array{Float64} = zeros(size(Slip[1,:]))

    taubefore::Array{Float64,2} = zeros(size(Slip))
    tauafter::Array{Float64,2} = zeros(size(Slip))
    
    hypo::Array{Float64} =  zeros(size(Slip[1,:]))   # Hypocenter
    vhypo::Array{Float64} = zeros(size(Slip[1,:]))   # Velocity at hypocenter

    Vthres = 0.001 # event threshold
    slipstart = 0
    it = 1; it2 = 1
    delfref = zeros(size(Slip[:,1]))

    for i = 1:length(Slip[1,:])

        # Start of each event
        if Vfmax[i] > 1.01*Vthres && slipstart == 0
            delfref = Slip[:,i]
            slipstart = 1
            tStart[it2] = time_[i]
            
            taubefore[:,it2] = Stress[:,i]
            vhypo[it2], indx = findmax(SlipVel[:,i])

            hypo[it2] = S.FltX[indx]

            it2 = it2+1
        end

        # End of each event
        if Vfmax[i] < 0.99*Vthres && slipstart == 1
            delfafter[:,it] = Slip[:,i] - delfref
            tauafter[:,it] = Stress[:,i]
            tEnd[it] = time_[i]
            slipstart = 0
            it = it + 1
        end
    end

    return delfafter[:,1:it-1], (taubefore-tauafter)[:,1:it-1], tStart[1:it2-1], tEnd[1:it-1], vhypo[1:it2-1], hypo[1:it2-1]
end

#..........................................................
# Compute the moment magnitude:
#       Assumed the rupture area to be square; the rupture
#       dimension along depth is the same as the rupture
#       dimension perpendicular to the plane
#..........................................................
function moment_magnitude_new(mu, FltX, delfafter, stressdrops ,time_)
    # Final coseismic slip of each earthquake
    #  delfafter, stressdrops = Coslip(S, Slip, SlipVel, Stress, time_)
    FltNglob = length(FltX)

    iter = length(delfafter[1,:])
    seismic_moment = zeros(iter)
    rupture_len = zeros(iter)
    fault_slip = zeros(iter)
    temp_sigma = 0
    iter2 = 1 

    del_sigma = zeros(iter)
    
    dx = diff(FltX).*1e3

    for i = 1:iter
        
        # slip threshold = 1% of maximum slip
        slip_thres = 0.01*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0; temp_sigma = 0; temp_slip = 0

        for j = 1:FltNglob
            if delfafter[j,i] >= slip_thres
                area = area + delfafter[j,i]*dx[j-1]
                zdim = zdim + dx[j-1]
                temp_slip = temp_slip + delfafter[j,i]

                # Avg. stress drops along rupture area
                temp_sigma = temp_sigma + stressdrops[j,i]*dx[j-1]
            end
        end
        
        seismic_moment[i] = mu*area*zdim
        del_sigma[i] = temp_sigma/zdim
        fault_slip[i] = temp_slip/zdim

        rupture_len[i] = zdim


    end
    #  seismic_moment = filter!(x->x!=0, seismic_moment)
    #  del_sigma = filter!(x->x!=0, del_sigma)
    Mw = (2/3)*log10.(seismic_moment.*1e7) .- 10.7

    return Mw, del_sigma, fault_slip, rupture_len, seismic_moment
end

