using JLD2
#
folder = "/mnt/d/"
#@load folder * "results_workspace.jld2" dat
#@load folder * "results_workspace_P.jld2" P
#@load folder * "results_workspace_datorig.jld2" dat_orig
#Pv = P

#GET PRESENCE/ABSENCE OF TARGETS REACHED

#Get measures of impact
function get_impact(dat_pre_int, dat_post_int)
    if !(dat_pre_int==0.) #can't have a 0 in denominator
        1 - dat_post_int/dat_pre_int #will get %decrease and %increase here
    else
        1 - dat_post_int/(dat_pre_int + 0.0000001) 
    end
end

function map_impacts(dat_pre_col, dat_post_col)
    map(x -> get_impact(dat_pre_col[x], dat_post_col[x]), 1:length(dat_pre_col))
end

#Using RCall to get % reaching target
using RCall
R"source('Scripts/libraries.R')"

#Fig 1 B
#Bar plot of RH, RA, and RE values - deterministic results
#initial parameter values
#          LH   LA   gH     gA     bHH         bAA         bHA         bAH         bAE         bEA         bEH         bHE         muH  muA  muE
p_B =     [0.1, 0.1, 0.001, 0.001, 0.1,        0.1,        0.1,        0.1,        0.1,        0.01,       0.01,       0.1,        0.1, 0.1, 0.2]
p_Bd =    [0.1, 0.1, 0.001, 0.001, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.1, 0.1, 0.2]
p_E =     [0.1, 0.1, 0.001, 0.001, 0.001,      0.001,      0.001,      0.001,      0.1420501,  0.1420501,  0.1420501,  0.1420501,  0.1, 0.1, 0.2]
p_A =     [0.1, 0.1, 0.001, 0.001, 0.001,      0.2019663,  0.001,      0.2019663,  0.2019663,  0.001,      0.001,      0.001,      0.1, 0.1, 0.2]
p_H =     [0.1, 0.1, 0.001, 0.001, 0.2019663,  0.001,      0.2019663,  0.001,      0.001,      0.001,      0.001,      0.2019663,  0.1, 0.1, 0.2]
p_orig =  [0.1, 0.1, 0.00,  0.00,  0.1,        0.1,        0.1,        0.1,        0.,         0.0,        0.0,        0.,         0.1, 0.1, 0.2]
p = [p_B, p_Bd, p_E, p_A, p_H, p_orig]

#       LH   LA   gH     gA     bHH         bAA         bHA         bAH         bAE         bEA         bEH         bHE         muH  muA  muE
p_B2 =  [0.1, 0.1, 0.001, 0.001, 0.1,        0.1,        0.1,        0.1,        0.1,        0.01,       0.01,       0.1,        0.1, 0.1, 0.2]
p_Bd2 = [0.1, 0.1, 0.001, 0.001, 0.08109928, 0.08109928, 0.08109928, 0.08109928, 0.08109928, 0.08109928, 0.08109928, 0.08109928, 0.1, 0.1, 0.2]
p_E2 =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.001,      0.001,      0.001,      0.23084954, 0.23084954, 0.23084954, 0.23084954, 0.1, 0.1, 0.2]
p_A2 =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.20239149, 0.001,      0.20239149, 0.20239149,  0.001,     0.001,      0.001,      0.1, 0.1, 0.2]
p_H2 =  [0.1, 0.1, 0.001, 0.001, 0.20239149, 0.001,      0.20239149,  0.001,      0.001,      0.001,     0.001,      0.20239149, 0.1, 0.1, 0.2]
p2 = [p_B2, p_Bd2, p_E2, p_A2, p_H2]

p_intervention = deepcopy(p)
foreach(xs -> xs[2] = 0, p_intervention)

include("model.jl")
R = map(x -> solve(ODEProblem(unboundeds, u0, tspan, p[x]))(1000), 1:6)
R2 = map(x -> solve(ODEProblem(boundeds, u0, tspan, p2[x]))(1000), 1:5)

R_int = map(x -> solve(ODEProblem(unboundeds, u0, tspan, p_intervention[x]))(1000), 1:6)

impact = map(x -> 1 - (R_int[x][1] / R[x][1]), 1:6)

R = R[1:5]
@rput R R2
R"source('Scripts/Fig 1/RHforeachTSplot.R')"

#Fig 2 A and B
#Compare bEH and LA interventions

#Transmission scenarios of interest: all
@load folder * "results_workspace_B_1.jld2" dat_B
@load folder * "results_workspace_Bd_1.jld2" dat_Bd
@load folder * "results_workspace_H_1.jld2" dat_H
@load folder * "results_workspace_A_1.jld2" dat_A
@load folder * "results_workspace_E_1.jld2" dat_E
dat_1 = [dat_B, dat_Bd, dat_H, dat_A, dat_E]

@load folder * "results_workspace_B_2.jld2" dat_B
@load folder * "results_workspace_Bd_2.jld2" dat_Bd
@load folder * "results_workspace_H_2.jld2" dat_H
@load folder * "results_workspace_A_2.jld2" dat_A
@load folder * "results_workspace_E_2.jld2" dat_E
dat_2 = [dat_B, dat_Bd, dat_H, dat_A, dat_E]

dat_B = Nothing
dat_Bd = Nothing
dat_H = Nothing
dat_A = Nothing
dat_E = Nothing

#Experiments needed: 
#bEH fixed to 0.1, LA fixed to 0.1 (exp. 4)
#bEH fixed to 0, LA fixed to 0.1 (exp. 6)
#bEH at TS, LA fixed to 0.1 (exp. 1)
#bEH at TS, LA fixed to 0 (exp. 2)
bEH_impact_unbounded = [map_impacts(dat_1[i][4][:,2], dat_1[i][6][:,2]) for i in 1:5]
bEH_impact_bounded = [map_impacts(dat_2[i][4][:,2], dat_2[i][6][:,2]) for i in 1:5]

LA_impact_unbounded = [map_impacts(dat_1[i][1][:,2], dat_1[i][2][:,2]) for i in 1:5]
LA_impact_bounded = [map_impacts(dat_2[i][1][:,2], dat_2[i][2][:,2]) for i in 1:5]

@rput bEH_impact_unbounded bEH_impact_bounded LA_impact_unbounded LA_impact_bounded
R"source('Scripts/Fig 2/fig2.R')"

#Fig 3 A, B and apendix
#Heatmaps of impact as LA and bEH changes, for two levels of bHA

#For appendix - all transmission scenarios, bounded and unbounded model, normal value of bHA
#Experiments needed
#bEH varying (unif), LA varying (unif), regular bHA - exp. 10
#bEH varying (unif), LA 0, regular bHA - exp. 9

heatmap_impacts_unbounded = [map_impacts(dat_1[i][10][:,2], dat_1[i][9][:,2]) for i in 1:5]
heatmap_impacts_bounded = [map_impacts(dat_2[i][10][:,2], dat_2[i][9][:,2]) for i in 1:5]

#Now for the Fig 3A and B - 
#Experiments needed
#bEH varying (unif), LA varying (unif), low bHA - exp. 8
#bEH varying, LA 0 - exp. 11, low bHA - exp. 7
#Only looking at balanced results

fig3A_impacts_unbounded = map_impacts(dat_1[4][8][:,2], dat_1[4][7][:,2])
fig3A_impacts_bounded = map_impacts(dat_2[4][8][:,2], dat_2[4][7][:,2])

include("parameters.jl")

@rput heatmap_impacts_unbounded heatmap_impacts_bounded fig3A_impacts_unbounded fig3A_impacts_bounded bEH_unif LA_unif
R"source('Scripts/Fig 3/Fig3A.R')"

#Fig 3C
#dotplot and lineplot of omegaA as betaEH increases
#Experiments needed
#varying bEH, LA fixed to 0.1 - exp 13
#varying bEH, LA fixed to 0 - exp.9 
fig3C_impacts_unbounded = [map_impacts(dat_1[i][13][:,2], dat_1[i][9][:,2]) for i in 1:5]
fig3C_impacts_bounded = [map_impacts(dat_2[i][13][:,2], dat_2[i][9][:,2]) for i in 1:5]

#Also need the results of the original model
@load folder * "results_workspace_datorig.jld2" dat_orig_unbounded
@load folder * "results_workspace_datorig_bo.jld2" dat_orig_bounded

impact_orig = map_impacts(dat_orig_unbounded[1][:,2], dat_orig_unbounded[2][:,2])
impact_orig_bounded = map_impacts(dat_orig_bounded[1][:,2], dat_orig_bounded[2][:,2])

@rput fig3C_impacts_bounded fig3C_impacts_unbounded impact_orig impact_orig_bounded
R"source('Scripts/Fig 3/Fig3C.R')"

#Fig 3 alternatively
#Point and error bar plot of omegaA in each transmission scenarios
#Fixed bEH (to each transmission scenario), LA fixed to 0.1 - experiment 1
#As above but LA fixed to 0 - experiment 2

f3_impacts_unbounded = [map_impacts(dat_1[i][1][:,2], dat_1[i][2]) for i in 1:5]
f3_impacts_bounded = [map_impacts(dat_2[i][1][:,2], dat_2[i][2]) for i in 1:5]
@rput f3_impacts_unbounded f3_impacts_bounded
R"source('Scripts/Fig 3/Fig3_alt.R')"
