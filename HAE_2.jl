using DifferentialEquations
using Distributions
using Plots
using RCall

# model function
function unboundeds(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, βHH, βAA, βHA, βAH, βAE, βEA, βEH, βHE, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + βHH*RH + βAH*RA + βEH*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βAA*RA + βHA*RH + βEA*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βAE*RA + βHE*RH - μE*RE
end

#check is working - us the same as the mathematica output to 6 decimal places at least
u0 = [0.0; 0.0; 0.0]
tspan = (0.0, 1000.)
p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.001, 0.1, 0.1, 0.01, 0.01, 0.1, 0.1, 0.1, 0.2]
prob = ODEProblem(unboundeds, u0, tspan, p)
sol = solve(prob)
plotlyjs() #plotting backend
plot(sol, ylims=(0.,1.), yticks=0.:.1:1.)

#sample paramter space - for all transmission scenarios.
#Make a little library for the distributions and parameters for the 4 varying parameters
#assume will always use lognormal for now...
struct lnp
    σ::Float64
    μ::Float64
end

struct bp
    α::Float64
    β::Float64
end

struct ts
    B
    H
    E
    A
end

struct vp
    βEH
    ΛH
    μE
    μH
end

#σ=1.
#pv = vp(ts(lnp(σ, 0.01), lnp(σ, 0.001),
#            lnp(σ, 0.14), lnp(σ, 0.001)), #βEH
#          lnp(σ, 0.1), #ΛH
#          lnp(σ, 0.2), #μE
#          lnp(σ, 0.1)) #μH

#beta distribution with mean = parameter selected in transmission scenarios, and var = mean/20
#parameters estimated from R code
pv = vp(ts(bp(0.188, 18.612), bp(0.01898, 18.96102),
          bp(2.268, 13.932), bp(0.01898, 18.96102)), #βEH
        bp(1.7, 15.3), #ΛH
        bp(3, 12), #μE
        bp(1.7, 15.3)) #μH


#non varying parameters
struct f_p
    ΛA::Float64
    γH::Float64
    γA::Float64
    βHH::ts
    βAA::ts
    βHA::ts
    βAH::ts
    βAE::ts
    βEA::ts
    βHE::ts
    μA::Float64
end

pf = f_p(0.1,0.001,0.001, #ΛA, γH, γA
         ts(0.1,0.2,0.001,0.001), #βHH
         ts(0.1,0.001,0.001,0.2), #βAA
         ts(0.001,0.2,0.001,0.001), #βHA
         ts(0.1,0.001,0.001,0.2), #βAH
         ts(0.1, 0.001, 0.14, 0.2), #βAE
         ts(0.01, 0.001, 0.14, 0.001), #βEA
         ts(0.1, 0.2, 0.14, 0.001), #βHE
         0.1) #μA

#set up parameter set that never changes first
N = 2000000
p_initial = zeros(N, 15)

#Parameters that have sampling distributions
#p_initial[:,1] .= rand(LogNormal(log(pv.ΛH.μ)+pv.ΛH.σ, pv.ΛH.σ), N)
#p_initial[:,15] .= rand(LogNormal(log(pv.μE.μ)+pv.μE.σ, pv.μE.σ), N)
#p_initial[:,13] .= rand(LogNormal(log(pv.μH.μ)+pv.μH.σ, pv.μH.σ), N)
p_initial[:,1] .= rand(Beta(pv.ΛH.α, pv.ΛH.β), N)
p_initial[:,15] .= rand(Beta(pv.μE.α, pv.μE.β), N)
p_initial[:,13] .= rand(Beta(pv.μH.α, pv.μH.β), N)

#Parameters that are fixed
#ΛA, γH, γA, βHH, βAA, βHA, βAH, βAE, βEA, βHE, μA
p_initial[:,[2,3,4,14]] .= [pf.ΛA, pf.γH, pf.γA, pf.μA]'

#Make parameter sets for the different transmission scenarios
p_initial_B = copy(p_initial)
p_initial_H = copy(p_initial)
p_initial_A = copy(p_initial)
p_initial_E = copy(p_initial)

#Parameters that have sampling distributions
#p_initial_B[:,11] .= rand(LogNormal(log(pv.βEH.B.μ)+pv.βEH.B.σ, pv.βEH.B.σ), N)
#p_initial_H[:,11] .= rand(LogNormal(log(pv.βEH.H.μ)+pv.βEH.H.σ, pv.βEH.H.σ), N)
#p_initial_A[:,11] .= rand(LogNormal(log(pv.βEH.A.μ)+pv.βEH.A.σ, pv.βEH.A.σ), N)
#p_initial_E[:,11] .= rand(LogNormal(log(pv.βEH.E.μ)+pv.βEH.E.σ, pv.βEH.E.σ), N)
p_initial_B[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_H[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_A[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_E[:,11] .= rand(Uniform(0.000001, 1.), N)

#Parameters that are fixed
#βHH, βAA, βHA, βAH, βAE, βEA, βHE
p_initial_B[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).B for x in [4,5,6,7,8,9,10]]'
p_initial_H[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).H for x in [4,5,6,7,8,9,10]]'
p_initial_A[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).A for x in [4,5,6,7,8,9,10]]'
p_initial_E[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).E for x in [4,5,6,7,8,9,10]]'

#get rid of negative numbers, or values over 1.5
function keep_ps(p)
    sets_LH = p[:,1] .< 1.
    sets_muH = p[:,13] .< 1.
    sets_muE = p[:,15] .<1.5
    keep = findall(sets_LH .& sets_muH .& sets_muE)
    return p[keep,:]
end
p_B = keep_ps(p_initial_B)
p_H = keep_ps(p_initial_H)
p_A = keep_ps(p_initial_A)
p_E = keep_ps(p_initial_E)

plot(
    histogram(p_B[:,[1,11,13,15]],
              xticks = range(0, 1.; step =0.1),xlims = (0.0, 1.),label = ["LH" "bEH" "mH" "mE"]),
    histogram(p_H[:,[1,11,13,15]],
            xticks = range(0, 1.; step =0.1),xlims = (0.0, 1.),label = ["LH" "bEH" "mH" "mE"]),
    histogram(p_A[:,[1,11,13,15]],
              xticks = range(0, 1; step =0.1),xlims = (0.0, 1.),label = ["LH" "bEH" "mH" "mE"]),
    histogram(p_E[:,[1,11,13,15]],
            xticks = range(0, 1.; step =0.1),xlims = (0.0, 1.),label = ["LH" "bEH" "mH" "mE"]),
    layout = 4)

#now numerically solve for each parameter set
#function for running model
function model_run(p, mod)
    dat = zeros(size(p)[1],2)
    u0 = [0.0;0.0;0.0]
    tspan = (0., 500.)

    #run 1
    @time for i in 1:size(p)[1]
      prob = ODEProblem(mod, u0, tspan, p[i,:])
      sol = solve(prob)
      dat[i, :] .= Array{Float64,1}(map(n -> sol(n)[1], [400,500]))
    end

    #run 2 - rerun for those possibly not at equilibrium
    not_eqlm = [abs(dat[i,2] - dat[i,1]) > 0.0000001 for i in 1:size(p)[1]]
    rerun = findall(not_eqlm)
    tspan = (0.,10000.)
    @time for i in rerun
        prob = ODEProblem(mod, u0, tspan, p[i,:])
        sol = solve(prob)
        dat[i, :] .= Array{Float64,1}(map(n -> sol(n)[1], [1900,2000]))
    end
    return dat
end

#RUN MODELS
#Index for runs
# 1. = varying bEH, LA fixed to 0.1.
# 2. = varying bEH, LA fixed to 0.0.
# 3. = varying bEH, varying LA.
# 4. = fixed bEH, varying LA.
# 5. = fixed bEH = 0, fixed LA = 0.1
# 6. = fixed bEH, LA = 0.

dat_E_1 = model_run(p_E, unboundeds) #250.2,1.58
dat_B_1 = model_run(p_B, unboundeds) #256.87, 1.66
dat_H_1 = model_run(p_H, unboundeds) #248.08, 1.65
dat_A_1 = model_run(p_A, unboundeds) #245.40, 1.21

#2. varying bEH, fixed LA = 0.0
p_E_2 = copy(p_E)
p_B_2 = copy(p_B)
p_H_2 = copy(p_H)
p_A_2 = copy(p_A)
p_B_2[:,2] .= 0
p_H_2[:,2] .= 0
p_A_2[:,2] .= 0
p_E_2[:,2] .= 0

dat_E_2 = model_run(p_E_2, unboundeds) #262.72, 2.08
dat_B_2 = model_run(p_B_2, unboundeds) #339.46, 227.85
dat_H_2 = model_run(p_H_2, unboundeds) #248.46, 1.78
dat_A_2 = model_run(p_A_2, unboundeds) #282.00, 2.71

#3. varying LA, varying bEH
p_E_3 = copy(p_E)
p_B_3 = copy(p_B)
p_H_3 = copy(p_H)
p_A_3 = copy(p_A)
p_E_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_E)[1])
p_B_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_B)[1])
p_H_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_H)[1])
p_A_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_A)[1])

dat_E_3 = model_run(p_E_3, unboundeds) #249.09, 1.40
dat_B_3 = model_run(p_B_3, unboundeds) #612.57, 1.60
dat_H_3 = model_run(p_H_3, unboundeds) #248.13, 1.63
dat_A_3 = model_run(p_A_3, unboundeds) #249.70, 1.14

#4. fixed bEH, varying LA
p_E_4 = copy(p_E_3)
p_B_4 = copy(p_B_3)
p_H_4 = copy(p_H_3)
p_A_4 = copy(p_A_3)
p_E_4[:,11] .= 0.14
p_B_4[:,11] .= 0.01
p_H_4[:,11] .= 0.001
p_A_4[:,11] .= 0.001

dat_E_4 = model_run(p_E_4, unboundeds) #249.09, 1.40
dat_B_4 = model_run(p_B_4, unboundeds) #612.57, 1.60
dat_H_4 = model_run(p_H_4, unboundeds) #248.13, 1.63
dat_A_4 = model_run(p_A_4, unboundeds) #249.70, 1.14

#5. fixed bEH = 0, fixed LA = 0.1
p_E_5 = copy(p_E)
p_B_5 = copy(p_B)
p_H_5 = copy(p_H)
p_A_5 = copy(p_A)
p_E_5[:,11] .= 0
p_B_5[:,11] .= 0
p_H_5[:,11] .= 0
p_A_5[:,11] .= 0
dat_E_5 = model_run(p_E_5, unboundeds) #262.72, 2.08
dat_B_5 = model_run(p_B_5, unboundeds) #339.46, 227.85
dat_H_5 = model_run(p_H_5, unboundeds) #248.46, 1.78
dat_A_5 = model_run(p_A_5, unboundeds) #282.00, 2.71

#6. fixed bEH, fixed LA = 0
p_E_6 = copy(p_E_4)
p_B_6 = copy(p_B_4)
p_H_6 = copy(p_H_4)
p_A_6 = copy(p_A_4)
p_B_6[:,2] .= 0
p_H_6[:,2] .= 0
p_A_6[:,2] .= 0
p_E_6[:,2] .= 0
dat_E_6 = model_run(p_E_6, unboundeds) #262.72, 2.08
dat_B_6 = model_run(p_B_6, unboundeds) #339.46, 227.85
dat_H_6 = model_run(p_H_6, unboundeds) #248.46, 1.78
dat_A_6 = model_run(p_A_6, unboundeds) #282.00, 2.71

#GET PRESENCE/ABSENCE OF TARGETS REACHED

#Get measures of impact
#1. for fixed LA 0.1 -> 0.0
function get_impact(dat1, dat2)
    if !(dat1==0.) #can't have a 0 in denominator
        1 - dat2/dat1 #will get %decrease and %increase here
    else
        1 - dat2/(dat1 + 0.0000001) #will get %decrease and %increase here
    end
end

#1.fixed LA intervention, varying bEH
impact_E = [get_impact(dat_E_1[i,2], dat_E_2[i,2]) for i in 1:size(p_E)[1]]
impact_B = [get_impact(dat_B_1[i,2], dat_B_2[i,2]) for i in 1:size(p_B)[1]]
impact_H = [get_impact(dat_H_1[i,2], dat_H_2[i,2]) for i in 1:size(p_H)[1]]
impact_A = [get_impact(dat_A_1[i,2], dat_A_2[i,2]) for i in 1:size(p_A)[1]]

#2. varyng LA intervention, varying bEH
impact_E_2 = [get_impact(dat_E_3[i,2], dat_E_2[i,2]) for i in 1:size(p_E)[1]]
impact_B_2 = [get_impact(dat_B_3[i,2], dat_B_2[i,2]) for i in 1:size(p_B)[1]]
impact_H_2 = [get_impact(dat_H_3[i,2], dat_H_2[i,2]) for i in 1:size(p_H)[1]]
impact_A_2 = [get_impact(dat_A_3[i,2], dat_A_2[i,2]) for i in 1:size(p_A)[1]]

#3. varying bEH intervention, fixed LA
impact_bEH_E = [get_impact(dat_E_1[i,2], dat_E_5[i,2]) for i in 1:size(p_E)[1]]
impact_bEH_B = [get_impact(dat_B_1[i,2], dat_B_5[i,2]) for i in 1:size(p_B)[1]]
impact_bEH_A = [get_impact(dat_A_1[i,2], dat_A_5[i,2]) for i in 1:size(p_A)[1]]
impact_bEH_H = [get_impact(dat_H_1[i,2], dat_H_5[i,2]) for i in 1:size(p_H)[1]]

#4. varying LA intervention, fixed bEH
impact_E_3 = [get_impact(dat_E_4[i,2], dat_E_6[i,2]) for i in 1:size(p_E)[1]]
impact_B_3 = [get_impact(dat_B_4[i,2], dat_B_6[i,2]) for i in 1:size(p_B)[1]]
impact_H_3 = [get_impact(dat_H_4[i,2], dat_H_6[i,2]) for i in 1:size(p_H)[1]]
impact_A_3 = [get_impact(dat_A_4[i,2], dat_A_6[i,2]) for i in 1:size(p_A)[1]]

#Did simulations reach target RH of 0.65 - 0.75?
n_target_B = [ifelse(0.65 < dat_B[i,2] <0.75, 1, 0) for i in 1:size(p_B)[1]]
n_target_H = [ifelse(0.65 < dat_H[i,2] <0.75, 1, 0) for i in 1:size(p_H)[1]]
n_target_A = [ifelse(0.65 < dat_A[i,2] <0.75, 1, 0) for i in 1:size(p_A)[1]]
n_target_E = [ifelse(0.65 < dat_E[i,2] <0.75, 1, 0) for i in 1:size(p_E)[1]]

#Did simulations have a low impact of less than 2%?
n_lowimpact_E = [ifelse(0. < impact_E[i] < 0.02,1,0) for i in 1:size(p_E)[1]]
n_lowimpact_B = [ifelse(0. < impact_B[i] < 0.02,1,0) for i in 1:size(p_B)[1]]
n_lowimpact_H = [ifelse(0. < impact_H[i] < 0.02,1,0) for i in 1:size(p_H)[1]]
n_lowimpact_A = [ifelse(0. < impact_A[i] < 0.02,1,0) for i in 1:size(p_A)[1]]

#What % of simulations reached the targets of interest?
#Bins for parameters
lower_bin = [0.:0.05:1;]
bin_N = size(lower_bin)[1]

#Using RCall to get % reaching target
R"""
source("M:/Github/animal-human-env-model/RGetPercAndPlot.R")
"""

@rput lower_bin p_E p_A p_B p_H p_E_3 p_A_3 p_B_3 p_H_3 p_E_4 p_A_4 p_B_4 p_H_4
@rput n_target_A n_target_B n_target_E n_target_H n_lowimpact_A n_lowimpact_B n_lowimpact_E n_lowimpact_H
@rput impact_A impact_B impact_E impact_H impact_A_2 impact_B_2 impact_E_2 impact_H_2 impact_A_3 impact_B_3 impact_E_3 impact_H_3
@rput impact_bEH_A impact_bEH_B impact_bEH_E impact_bEH_H

# Conclusion 1: realistic RHs are attainable for environmental transmission scenarios
R"""
bEHmuE_df_E = get_perc_target(lower_bin, p_E[, c(11, 15)], n_target_E)
bEHmuH_df_E = get_perc_target(lower_bin, p_E[, c(11, 13)], n_target_E)
bEHLH_df_E = get_perc_target(lower_bin, p_E[, c(11, 1)], n_target_E)

bEHmuE_df_B = get_perc_target(lower_bin, p_B[, c(11, 15)], n_target_B)
bEHmuH_df_B = get_perc_target(lower_bin, p_B[, c(11, 13)], n_target_B)
bEHLH_df_B = get_perc_target(lower_bin, p_B[, c(11, 1)], n_target_B)

bEHmuE_df_A = get_perc_target(lower_bin, p_A[, c(11, 15)], n_target_A)
bEHmuH_df_A = get_perc_target(lower_bin, p_A[, c(11, 13)], n_target_A)
bEHLH_df_A = get_perc_target(lower_bin, p_A[, c(11, 1)], n_target_A)

bEHmuE_df_H = get_perc_target(lower_bin, p_H[, c(11, 15)], n_target_H)
bEHmuH_df_H = get_perc_target(lower_bin, p_H[, c(11, 13)], n_target_H)
bEHLH_df_H = get_perc_target(lower_bin, p_H[, c(11, 1)], n_target_H)
"""

#plotting!
R"""
p1 = plot_heatmap(bEHmuE_df_E, c('bEH', 'muE'), '% Target Achieved, ts = E', '', limits = c(0, 0.75))
p2 = plot_heatmap(bEHmuH_df_E, c('bEH', 'muH'), '% Target Achieved, ts = E', '', limits = c(0, 0.75))
p3 = plot_heatmap(bEHLH_df_E, c('bEH', 'LH'), '% Target Achieved, ts = E', '', limits = c(0, 0.75))

p4 = plot_heatmap(bEHmuE_df_B, c('bEH', 'muE'), '% Target Achieved, ts = B', '', limits = c(0, 0.75))
p5 = plot_heatmap(bEHmuH_df_B, c('bEH', 'muH'), '% Target Achieved, ts = B', '', limits = c(0, 0.75))
p6 = plot_heatmap(bEHLH_df_B, c('bEH', 'LH'), '% Target Achieved, ts = B', '', limits = c(0, 0.75))

p7 = plot_heatmap(bEHmuE_df_H, c('bEH', 'muE'), '% Target Achieved, ts = H', '', limits = c(0, 0.75))
p8 = plot_heatmap(bEHmuH_df_H, c('bEH', 'muH'), '% Target Achieved, ts = H', '', limits = c(0, 0.75))
p9 = plot_heatmap(bEHLH_df_H, c('bEH', 'LH'), '% Target Achieved, ts = H', '', limits = c(0, 0.75))

p10 = plot_heatmap(bEHmuE_df_A, c('bEH', 'muE'), '% Target Achieved, ts = A', '', limits = c(0, 0.75))
p11 = plot_heatmap(bEHmuH_df_A, c('bEH', 'muH'), '% Target Achieved, ts = A', '', limits = c(0, 0.75))
p12 = plot_heatmap(bEHLH_df_A, c('bEH', 'LH'), '% Target Achieved, ts = A', '', limits = c(0, 0.75))
"""

R"""
#svg('M:/Project folders/Model env compartment/Plots/ptaplot.svg', height=12, width = 20)
ptaplot <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, nrow = 4, ncol = 3)
#dev.off()
"""

#Conclusions 2: impact of reducing LA is low for parameter combinations of interest
R"""
bEHmuE_df_E_low_impact = get_perc_target(lower_bin, p_E[, c(11, 15)], n_lowimpact_E)
bEHmuH_df_E_low_impact = get_perc_target(lower_bin, p_E[, c(11, 13)], n_lowimpact_E)
bEHLH_df_E_low_impact = get_perc_target(lower_bin, p_E[, c(11, 1)], n_lowimpact_E)

bEHmuE_df_B_low_impact = get_perc_target(lower_bin, p_B[, c(11, 15)], n_lowimpact_B)
bEHmuH_df_B_low_impact = get_perc_target(lower_bin, p_B[, c(11, 13)], n_lowimpact_B)
bEHLH_df_B_low_impact = get_perc_target(lower_bin, p_B[, c(11, 1)], n_lowimpact_B)

bEHmuE_df_A_low_impact = get_perc_target(lower_bin, p_A[, c(11, 15)], n_lowimpact_A)
bEHmuH_df_A_low_impact = get_perc_target(lower_bin, p_A[, c(11, 13)], n_lowimpact_A)
bEHLH_df_A_low_impact = get_perc_target(lower_bin, p_A[, c(11, 1)], n_lowimpact_A)

bEHmuE_df_H_low_impact = get_perc_target(lower_bin, p_H[, c(11, 15)], n_lowimpact_H)
bEHmuH_df_H_low_impact = get_perc_target(lower_bin, p_H[, c(11, 13)], n_lowimpact_H)
bEHLH_df_H_low_impact = get_perc_target(lower_bin, p_H[, c(11, 1)], n_lowimpact_H)
"""

R"""
p13 = plot_heatmap(bEHmuE_df_E_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = E', '', limits = c(0,1))
p14 = plot_heatmap(bEHmuH_df_E_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = E', '', limits = c(0,1))
p15 = plot_heatmap(bEHLH_df_E_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = E', '', limits = c(0,1))

p16 = plot_heatmap(bEHmuE_df_B_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = B', '', limits = c(0,1))
p17 = plot_heatmap(bEHmuH_df_B_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = B', '', limits = c(0,1))
p18 = plot_heatmap(bEHLH_df_B_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = B', '', limits = c(0,1))

p19 = plot_heatmap(bEHmuE_df_H_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = H', '', limits = c(0,1))
p20 = plot_heatmap(bEHmuH_df_H_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = H', '', limits = c(0,1))
p21 = plot_heatmap(bEHLH_df_H_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = H', '', limits = c(0,1))

p22 = plot_heatmap(bEHmuE_df_A_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = A', '', limits = c(0,1))
p23 = plot_heatmap(bEHmuH_df_A_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = A', '', limits = c(0,1))
p24 = plot_heatmap(bEHLH_df_A_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = A', '', limits = c(0,1))
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/liplot.svg', height=12, width = 20)
liplot = grid.arrange(p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,nrow = 4, ncol = 3)
#dev.off()
"""


#Conclusion 3: increasing βEH in many cases reduces the impacts of increased LA
R"""
bEHmuE_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 15)], impact_E)
bEHmuH_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 13)], impact_E)
bEHLH_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 1)], impact_E)

bEHmuE_df_B_impact = get_mean_var(lower_bin, p_B[, c(11, 15)], impact_B)
bEHmuH_df_B_impact = get_mean_var(lower_bin, p_B[, c(11, 13)], impact_B)
bEHLH_df_B_impact = get_mean_var(lower_bin, p_B[, c(11, 1)], impact_B)

bEHmuE_df_A_impact = get_mean_var(lower_bin, p_A[, c(11, 15)], impact_A)
bEHmuH_df_A_impact = get_mean_var(lower_bin, p_A[, c(11, 13)], impact_A)
bEHLH_df_A_impact = get_mean_var(lower_bin, p_A[, c(11, 1)], impact_A)

bEHmuE_df_H_impact = get_mean_var(lower_bin, p_H[, c(11, 15)], impact_H)
bEHmuH_df_H_impact = get_mean_var(lower_bin, p_H[, c(11, 13)], impact_H)
bEHLH_df_H_impact = get_mean_var(lower_bin, p_H[, c(11, 1)], impact_H)
"""

R"""
p25_1 = plot_heatmap(bEHmuE_df_E_impact, c('bEH', 'muE'), 'Mean impact, ts = E', '', limits = c(0, 0.4))
p25 = plot_heatmap(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Mean impact, ts = E', '', limits = c(0, 0.4))
p26 = plot_heatmap(bEHLH_df_E_impact, c('bEH', 'LH'), 'Mean impact, ts = E', '', limits = c(0, 0.4))

p27 = plot_heatmap(bEHmuE_df_B_impact, c('bEH', 'muE'), 'Mean impact, ts = B', '', limits = c(0, 0.4))
p28 = plot_heatmap(bEHmuH_df_B_impact, c('bEH', 'muH'), 'Mean impact, ts = B', '', limits = c(0, 0.4))
p29 = plot_heatmap(bEHLH_df_B_impact, c('bEH', 'LH'), 'Mean impact, ts = B', '', limits = c(0, 0.4))

p30 = plot_heatmap(bEHmuE_df_H_impact, c('bEH', 'muE'), 'Mean impact, ts = H', '', limits = c(0, 0.4))
p31 = plot_heatmap(bEHmuH_df_H_impact, c('bEH', 'muH'), 'Mean impact, ts = H', '', limits = c(0, 0.4))
p32 = plot_heatmap(bEHLH_df_H_impact, c('bEH', 'LH'), 'Mean impact, ts = H', '', limits = c(0, 0.4))

p33 = plot_heatmap(bEHmuE_df_A_impact, c('bEH', 'muE'), 'Mean impact, ts = A', '', limits = c(0, 0.4))
p34 = plot_heatmap(bEHmuH_df_A_impact, c('bEH', 'muH'), 'Mean impact, ts = A', '', limits = c(0, 0.4))
p35 = plot_heatmap(bEHLH_df_A_impact, c('bEH', 'LH'), 'Mean impact, ts = A', '', limits = c(0, 0.4))
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miplot.svg', height=12, width = 20)
miplot = grid.arrange(p25_1,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,nrow=4,ncol=3)
#dev.off()
"""

R"""
p36 = plot_heatmap_var(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Impact variance, ts = E', '', limits = c(0, 0.1))
p37 = plot_heatmap_var(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Impact variance, ts = E', '', limits = c(0, 0.1))
p38 = plot_heatmap_var(bEHLH_df_E_impact, c('bEH', 'LH'), 'Impact variance, ts = E', '', limits = c(0, 0.1))

p39 = plot_heatmap_var(bEHmuE_df_B_impact, c('bEH', 'muE'), 'Impact variance, ts = B', '', limits = c(0, 0.1))
p40 = plot_heatmap_var(bEHmuH_df_B_impact, c('bEH', 'muH'), 'Impact variance, ts = B', '', limits = c(0, 0.1))
p41 = plot_heatmap_var(bEHLH_df_B_impact, c('bEH', 'LH'), 'Impact variance, ts = B', '', limits = c(0, 0.1))

p42 = plot_heatmap_var(bEHmuE_df_H_impact, c('bEH', 'muE'), 'Impact variance, ts = H', '', limits = c(0, 0.1))
p43 = plot_heatmap_var(bEHmuH_df_H_impact, c('bEH', 'muH'), 'Impact variance, ts = H', '', limits = c(0, 0.1))
p44 = plot_heatmap_var(bEHLH_df_H_impact, c('bEH', 'LH'), 'Impact variance, ts = H', '', limits = c(0, 0.1))

p45 = plot_heatmap_var(bEHmuE_df_A_impact, c('bEH', 'muE'), 'Impact variance, ts = A', '', limits = c(0, 0.1))
p46 = plot_heatmap_var(bEHmuH_df_A_impact, c('bEH', 'muH'), 'Impact variance, ts = A', '', limits = c(0, 0.1))
p47 = plot_heatmap_var(bEHLH_df_A_impact, c('bEH', 'LH'), 'Impact variance, ts = A', '', limits = c(0, 0.1))

#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/viplot.svg', height=12, width = 20)
viplot = grid.arrange(p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,nrow=4,ncol=3)
#dev.off()
"""

R"""
impact_bEHLA_E = get_mean_var(lower_bin, p_E_3[, c(11, 2)], impact_E_2)
impact_bEHLA_H = get_mean_var(lower_bin, p_H_3[, c(11, 2)], impact_H_2)
impact_bEHLA_A = get_mean_var(lower_bin, p_A_3[, c(11, 2)], impact_A_2)
impact_bEHLA_B = get_mean_var(lower_bin, p_B_3[, c(11, 2)], impact_B_2)

p48 = plot_heatmap(impact_bEHLA_E, c('bEH', 'LA'), 'Mean impact, ts = E', '', limits = c(0.,max(impact_bEHLA_E$mean)))
p49 = plot_heatmap(impact_bEHLA_B, c('bEH', 'LA'), 'Mean impact, ts = B', '', limits = c(0.,max(impact_bEHLA_B$mean)))
p50 = plot_heatmap(impact_bEHLA_H, c('bEH', 'LA'), 'Mean impact, ts = H', '', limits = c(0.,max(impact_bEHLA_H$mean)))
p51 = plot_heatmap(impact_bEHLA_A, c('bEH', 'LA'), 'Mean impact, ts = A', '', limits = c(0.,max(impact_bEHLA_A$mean)))

p52 = plot_heatmap_var(impact_bEHLA_E, c('bEH', 'LA'), 'Impact variance, ts = E', '', limits = c(0.,max(impact_bEHLA_E$var)))
p53 = plot_heatmap_var(impact_bEHLA_B, c('bEH', 'LA'), 'Impact variance, ts = B', '', limits = c(0.,max(impact_bEHLA_B$var)))
p54 = plot_heatmap_var(impact_bEHLA_H, c('bEH', 'LA'), 'Impact variance, ts = H', '', limits = c(0.,max(impact_bEHLA_H$var)))
p55 = plot_heatmap_var(impact_bEHLA_A, c('bEH', 'LA'), 'Impact variance, ts = A', '', limits = c(0.,max(impact_bEHLA_A$var)))
"""

R"""
svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miLAplot.svg', height=12, width = 7)
miLAplot = grid.arrange(p48, p49, p50, p51, nrow=4)
dev.off()

#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/viLAplot.svg', height=12, width = 7)
#viLAplot = grid.arrange(p52, p53, p54, p55, nrow=4)
#dev.off()
"""

R"""
df_E <- data.frame(param_val = c(p_E[1:10000,11],  p_E_4[1:10000,2]),
                   param = rep(c("beH", "LA"), each = 10000),
                   impact = c(impact_bEH_E[1:10000], impact_E_3[1:10000]))
df_A <- data.frame(param_val = c(p_A[1:10000,11],  p_A_4[1:10000,2]),
                  param = rep(c("beH", "LA"), each = 10000),
                  impact = c(impact_bEH_A[1:10000], impact_A_3[1:10000]))
df_B <- data.frame(param_val = c(p_B[1:10000,11],  p_B_4[1:10000,2]),
                 param = rep(c("beH", "LA"), each = 10000),
                 impact = c(impact_bEH_B[1:10000], impact_B_3[1:10000]))
df_H <- data.frame(param_val = c(p_H[1:10000,11],  p_H_4[1:10000,2]),
                   param = rep(c("beH", "LA"), each = 10000),
                   impact = c(impact_bEH_H[1:10000], impact_H_3[1:10000]))

nice_plot <- function(df, ts) {
    ggplot(df, aes(param_val, impact, col = param)) +
                    geom_point(shape = 1) +
                    labs(x = "Parameter value",
                         y = "Impact",
                         col = "Parameter targeted \nby intervention",
                        title = paste("Transmission scenario: ", ts, sep = "")) +
                    geom_smooth(method = "loess", colour = "black") +
                    facet_grid(.~param) +
                    theme_bw()
}

p49 <- nice_plot(df_E, "Environment dominated")

p50 <- nice_plot(df_H, "Human dominated")

p51 <- nice_plot(df_A, "Animal dominated")

p52 <- nice_plot(df_B, "Balanced")

svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/bEHLAimpactcomparison.svg', height=12, width = 7)
bEHLAimpactcomparison = grid.arrange(p49, p50, p51, p52, nrow = 4)
dev.off()
"""

#What is the range of impacts?
#Get parameter sets with +- 10% of bEH value of that transmission scenario
bEH_E = findall((0.95 * 0.14) .< p_E[:,11] .< (1.05 * 0.14))
dat_E_bEH_target = dat_E[bEH_E,2]

bEH_B = findall((0.95 * 0.01) .< p_B[:,11] .< (1.05 * 0.01))
dat_B_bEH_target = dat_B[bEH_B,2]

bEH_A = findall((0.95 * 0.001) .< p_A[:,11] .< (1.05 * 0.001))
dat_A_bEH_target = dat_A[bEH_A,2]

bEH_H = findall((0.95 * 0.001) .< p_H[:,11] .< (1.05 * 0.001))
dat_H_bEH_target = dat_H[bEH_H,2]

#Get together RH measures for these indices
@rput dat_E_bEH_target dat_B_bEH_target dat_H_bEH_target dat_A_bEH_target

R"""
RH <- data.frame(ts = c(rep('E', length(dat_E_bEH_target)),
                        rep('B', length(dat_B_bEH_target)),
                        rep('A', length(dat_A_bEH_target)),
                        rep('H', length(dat_H_bEH_target))),
                 RH = c(dat_E_bEH_target,
                        dat_B_bEH_target,
                        dat_A_bEH_target,
                        dat_H_bEH_target))
ggplot(RH, aes(ts, RH)) + geom_boxplot() +
    geom_hline(yintercept = 0.71)
"""
