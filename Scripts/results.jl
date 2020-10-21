using JLD2
#
@load "D:/results_workspace.jld2" dat
@load "D:/results_workspace_P.jld2" P
@load "D:/results_workspace_datorig.jld2" dat_orig
Pv = P

#GET PRESENCE/ABSENCE OF TARGETS REACHED

#Get measures of impact
function get_impact(dat_pre_int, dat_post_int)
    if !(dat_pre_int==0.) #can't have a 0 in denominator
        1 - dat_post_int/dat_pre_int #will get %decrease and %increase here
    else
        1 - dat_post_int/(dat_pre_int + 0.0000001) #will get %decrease and %increase here
    end
end

function map_impacts(dat_pre_col, dat_post_col)
    map(x -> get_impact(dat_pre_col[x], dat_post_col[x]), 1:length(dat_pre_col))
end

#1.fixed LA intervention, varying bEH
impact_1 = [map_impacts(dat[i,1][:,2], dat[i,2][:,2]) for i in 1:5]

#2. varyng LA intervention, varying bEH
impact_2 = [map_impacts(dat[i,3][:,2], dat[i,2][:,2]) for i in 1:5]

#3. varying bEH intervention, fixed LA
impact_3 = [map_impacts(dat[i,1][:,2], dat[i,5][:,2]) for i in 1:5]

#4. varying LA intervention, fixed bEH
impact_4 = [map_impacts(dat[i,4][:,2], dat[i,6][:,2]) for i in 1:5]

#5. varying LA intervention, varying bEH, low bHA
impact_5 = [map_impacts(dat[i,9][:,2], dat[i, 8][:,2]) for i in 1:5]

#6. fixed bEH intervention
impact_6 = [map_impacts(dat[i,7][:,2], dat[i,5][:,2]) for i in 1:5]

#7. fixed LA intervention
impact_7 = [map_impacts(dat[i,7][:,2], dat[i,6][:,2]) for i in 1:5]

#8. both bEH and LA reduced to 0
impact_8 = [map_impacts(dat[i,7][:,2], dat[i, 10][:,2]) for i in 1:5]

#9. impact of original model - high impact scenario
impact_9 = map_impacts(dat_orig[1][:,2], dat_orig[2][:,2])

#Did simulations reach target RH of 0.65 - 0.75?
n_target = hcat([Int.(0.65 .< dat[x,1][:,1] .< 0.75) for x in 1:5])

#Did simulations have a low impact of less than 2%?
n_low_impact = hcat([Int.(impact_1[x] .< 0.02) for x in 1:5])


#Using RCall to get % reaching target
#rmprocs(11)
using RCall
R"source('Scripts/libraries.R')"

#Fig 1 B
#Bar plot of RH, RA, and RE values - deterministic results
#include("Scripts/model.jl")
#initial parameter values
#       LH   LA   gH     gA     bHH         bAA         bHA         bAH         bAE         bEA         bEH         bHE         muH  muA  muE
p_B =  [0.1, 0.1, 0.001, 0.001, 0.1,        0.1,        0.001,      0.1,        0.1,        0.01,       0.01,       0.1,        0.1, 0.1, 0.2]
p_Bd = [0.1, 0.1, 0.001, 0.001, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.1, 0.1, 0.2]
p_E =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.001,      0.001,      0.001,      0.1420501,  0.1420501,  0.1420501,  0.1420501,  0.1, 0.1, 0.2]
p_A =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.2019663,  0.001,      0.2019663,  0.2019663,  0.001,      0.001,      0.001,      0.1, 0.1, 0.2]
p_H =  [0.1, 0.1, 0.001, 0.001, 0.2019663,  0.001,      0.2019663,  0.001,      0.001,      0.001,      0.001,      0.2019663,  0.1, 0.1, 0.2]
p = [p_B, p_Bd, p_E, p_A, p_H]
R = map(x -> solve(ODEProblem(unboundeds, u0, tspan, p[x]))(200), 1:5)
@rput R
R"source('Scripts/Fig 1/RHforeachTSplot.R')"

#Fig 2 A and B 
#Compare bEH and LA interventions
Pv = P
@rput Pv impact_3 impact_4
R"source('Scripts/Fig 2/fig2.R')"

#Fig 3 A, B and apendix
#Heatmaps of impact as LA and bEH changes, for two levels of bHA
#Bins for parameters
lower_bin = [0.:0.05:1;]
bin_N = size(lower_bin)[1]
@rput impact_2 impact_5 lower_bin bin_N
R"source('Scripts/Fig 3/Fig3A.R')"

#Fig 3C
#dotplot and lineplot of omegaA as betaEH increases
@rput impact_1 impact_9
R"source('Scripts/Fig 3/Fig3C.R')"