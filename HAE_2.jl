using DifferentialEquations
using Distributions
using Plots
using CSV
using Tables

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

struct ts
    B
    H
    E
    A
end

struct vp
    βEH::ts
    ΛH::lnp
    μE::lnp
    μH::lnp
end

σ=1.
pv = vp(ts(lnp(σ, 0.01), lnp(σ, 0.001),
            lnp(σ, 0.14), lnp(σ, 0.001)), #βEH
          lnp(σ, 0.1), #ΛH
          lnp(σ, 0.2), #μE
          lnp(σ, 0.1)) #μH

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
p_initial[:,1] .= rand(LogNormal(log(pv.ΛH.μ)+pv.ΛH.σ, pv.ΛH.σ), N)
p_initial[:,15] .= rand(LogNormal(log(pv.μE.μ)+pv.μE.σ, pv.μE.σ), N)
p_initial[:,13] .= rand(LogNormal(log(pv.μH.μ)+pv.μH.σ, pv.μH.σ), N)
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
    keep=[]
    @time for i in 1:size(p)[1]
        if !(any(p[i,:] .< 0.))
            if !(any(p[i,:] .> 1.5))
                if !(p[i,13] > 1.) #no big values of μH
                    push!(keep, i)
                end
            end
        end
    end
    return p[keep, :]
end
p_B = keep_ps(p_initial_B)
p_H = keep_ps(p_initial_H)
p_A = keep_ps(p_initial_A)
p_E = keep_ps(p_initial_E)

plot(
    histogram(p_B[:,[1,11,13,15]],
              xticks = range(0, 1.5; step =0.1),xlims = (0.0, 1.5),label = ["ΛH" "βEH" "μH" "μE"]),
    histogram(p_H[:,[1,11,13,15]],
            xticks = range(0, 1.5; step =0.1),xlims = (0.0, 1.5),label = ["ΛH" "βEH" "μH" "μE"]),
    histogram(p_A[:,[1,11,13,15]],
              xticks = range(0, 1.5; step =0.1),xlims = (0.0, 1.5),label = ["ΛH" "βEH" "μH" "μE"]),
    histogram(p_E[:,[1,11,13,15]],
            xticks = range(0, 1.5; step =0.1),xlims = (0.0, 1.5),label = ["ΛH" "βEH" "μH" "μE"]),
    layout = 4)

#now numerically solve for each parameter set
#function for running model
function model_run(p, mod)
    dat = zeros(size(p)[1],2)
    u0 = [0.0;0.0;0.0]
    tspan = (0., 500.)
    re_run = []
    #run 1
    @time for i in 1:size(p)[1]
      prob = ODEProblem(mod, u0, tspan, p[i,:])
      sol = solve(prob)
      dat[i, :] .= Array{Float64,1}(map(n -> sol(n)[1], [400,500]))
      if abs(dat[i,2] - dat[i,1]) > 0.0000001
          push!(re_run, i)
      end
    end

    #rerun for those possibly not at equilibrium
    tspan = (0.,10000.)
    @time for i in re_run
        prob = ODEProblem(mod, u0, tspan, p[i,:])
        sol = solve(prob)
        dat[i, :] .= Array{Float64,1}(map(n -> sol(n)[1], [1900,2000]))
    end
    return dat
end
model_run(p_E[1:5,:], unboundeds) #short run to get the function going

dat_E = model_run(p_E, unboundeds)
dat_B = model_run(p_B, unboundeds)
dat_H = model_run(p_H, unboundeds)
dat_A = model_run(p_A, unboundeds)

#Find number of runs where 0.65 < RH < 0.75
n_target_B = [ifelse(0.65 < dat_B[i,2] <0.75, 1, 0) for i in 1:size(p_B)[1]]
n_target_H = [ifelse(0.65 < dat_H[i,2] <0.75, 1, 0) for i in 1:size(p_H)[1]]
n_target_A = [ifelse(0.65 < dat_A[i,2] <0.75, 1, 0) for i in 1:size(p_A)[1]]
n_target_E = [ifelse(0.65 < dat_E[i,2] <0.75, 1, 0) for i in 1:size(p_E)[1]]

#Bins for parameters
lower_bin = [0.:0.05:1;]
bin_N = size(lower_bin)[1]

#Using RCall to get % reaching target
R"""
source("M:/Project\ folders\\/Model\ env\ compartment\\/Scipts\\/RGetPercAndPlot.R")
"""

@rput lower_bin p_E p_A p_B p_H n_target_A n_target_B n_target_E n_target_H

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
p1 = plot_heatmap(bEHmuE_df_E, c('bEH', 'muE'), '% Target Achieved, ts = E', '')
p2 = plot_heatmap(bEHmuH_df_E, c('bEH', 'muH'), '% Target Achieved, ts = E', '')
p3 = plot_heatmap(bEHLH_df_E, c('bEH', 'LH'), '% Target Achieved, ts = E', '')

p4 = plot_heatmap(bEHmuE_df_B, c('bEH', 'muE'), '% Target Achieved, ts = B', '')
p5 = plot_heatmap(bEHmuH_df_B, c('bEH', 'muH'), '% Target Achieved, ts = B', '')
p6 = plot_heatmap(bEHLH_df_B, c('bEH', 'LH'), '% Target Achieved, ts = B', '')

p7 = plot_heatmap(bEHmuE_df_H, c('bEH', 'muE'), '% Target Achieved, ts = H', '')
p8 = plot_heatmap(bEHmuH_df_H, c('bEH', 'muH'), '% Target Achieved, ts = H', '')
p9 = plot_heatmap(bEHLH_df_H, c('bEH', 'LH'), '% Target Achieved, ts = H', '')

p10 = plot_heatmap(bEHmuE_df_A, c('bEH', 'muE'), '% Target Achieved, ts = A', '')
p11 = plot_heatmap(bEHmuH_df_A, c('bEH', 'muH'), '% Target Achieved, ts = A', '')
p12 = plot_heatmap(bEHLH_df_A, c('bEH', 'LH'), '% Target Achieved, ts = A', '')
"""

R"""
svg('M:/Project folders/Model env compartment/Plots/ptaplot.svg', height=12, width = 20)
ptaplot <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, nrow = 4, ncol = 3)
dev.off()
"""

#Conclusions 2: impact of reducing LA is low for parameter combinations of interest
#I think I should do this for human transmission scenario and then the animal transmission scenario - best and worst chance to have an impact

#run with 0 lambda value
p_B[:,2] .= 0
p_H[:,2] .= 0
p_A[:,2] .= 0
p_E[:,2] .= 0
dat_E_noLA = model_run(p_E, unboundeds)
dat_B_noLA = model_run(p_B, unboundeds)
dat_H_noLA = model_run(p_H, unboundeds)
dat_A_noLA = model_run(p_A, unboundeds)

impact_E = [ifelse(!(dat_E[i]==0.), #can't have a 0 in denominator, and if 0 wouldn't expect any impact anyway
                   1 - dat_E_noLA[i]/dat_E[i], #will get %decrease and %increase here, will maybe cut out later
                   0) for i in 1:size(p_E)[1]]
impact_B = [ifelse(!(dat_B[i]==0.),
                  1 - dat_B_noLA[i]/dat_B[i],
                  0) for i in 1:size(p_E)[1]]
impact_H = [ifelse(!(dat_H[i]==0.),
                 1 - dat_H_noLA[i]/dat_H[i],
                 0) for i in 1:size(p_E)[1]]
impact_A = [ifelse(!(dat_A[i]==0.),
                1 - dat_A_noLA[i]/dat_A[i],
                0) for i in 1:size(p_E)[1]]

n_lowimpact_E = [ifelse(0. < impact_E[i] < 0.02,
                        1,
                        0) for i in 1:size(p_E)[1]]
n_lowimpact_B = [ifelse(0. < impact_B[i] < 0.02,
                        1,
                        0) for i in 1:size(p_B)[1]]
n_lowimpact_H = [ifelse(0. < impact_H[i] < 0.02,
                        1,
                        0) for i in 1:size(p_H)[1]]
n_lowimpact_A = [ifelse(0. < impact_A[i] < 0.02,
                        1,
                        0) for i in 1:size(p_A)[1]]

@rput impact_E impact_B impact_A impact_H n_lowimpact_E n_lowimpact_B n_lowimpact_A n_lowimpact_H
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
p13 = plot_heatmap(bEHmuE_df_E_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = E', '')
p14 = plot_heatmap(bEHmuH_df_E_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = E', '')
p15 = plot_heatmap(bEHLH_df_E_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = E', '')

p16 = plot_heatmap(bEHmuE_df_B_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = B', '')
p17 = plot_heatmap(bEHmuH_df_B_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = B', '')
p18 = plot_heatmap(bEHLH_df_B_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = B', '')

p19 = plot_heatmap(bEHmuE_df_H_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = H', '')
p20 = plot_heatmap(bEHmuH_df_H_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = H', '')
p21 = plot_heatmap(bEHLH_df_H_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = H', '')

p22 = plot_heatmap(bEHmuE_df_A_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = A', '')
p23 = plot_heatmap(bEHmuH_df_A_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = A', '')
p24 = plot_heatmap(bEHLH_df_A_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = A', '')
"""

R"""
svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/liplot.svg', height=12, width = 20)
liplot = grid.arrange(p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,nrow = 4, ncol = 3)
dev.off()
"""


#Final conclusion - increasing βEH in many cases reduces the impacts of increased LH
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
p25_1 = plot_heatmap(bEHmuE_df_E_impact, c('bEH', 'muE'), 'Mean impact, ts = E', '')
p25 = plot_heatmap(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Mean impact, ts = E', '')
p26 = plot_heatmap(bEHLH_df_E_impact, c('bEH', 'LH'), 'Mean impact, ts = E', '')

p27 = plot_heatmap(bEHmuE_df_B_impact, c('bEH', 'muE'), 'Mean impact, ts = B', '')
p28 = plot_heatmap(bEHmuH_df_B_impact, c('bEH', 'muH'), 'Mean impact, ts = B', '')
p29 = plot_heatmap(bEHLH_df_B_impact, c('bEH', 'LH'), 'Mean impact, ts = B', '')

p30 = plot_heatmap(bEHmuE_df_H_impact, c('bEH', 'muE'), 'Mean impact, ts = H', '')
p31 = plot_heatmap(bEHmuH_df_H_impact, c('bEH', 'muH'), 'Mean impact, ts = H', '')
p32 = plot_heatmap(bEHLH_df_H_impact, c('bEH', 'LH'), 'Mean impact, ts = H', '')

p33 = plot_heatmap(bEHmuE_df_A_impact, c('bEH', 'muE'), 'Mean impact, ts = A', '')
p34 = plot_heatmap(bEHmuH_df_A_impact, c('bEH', 'muH'), 'Mean impact, ts = A', '')
p35 = plot_heatmap(bEHLH_df_A_impact, c('bEH', 'LH'), 'Mean impact, ts = A', '')
"""

R"""
svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miplot.svg', height=12, width = 20)
miplot = grid.arrange(p25_1,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,nrow=4,ncol=3)
dev.off()
"""

R"""
p36 = plot_heatmap_var(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Impact variance, ts = E', '')
p37 = plot_heatmap_var(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Impact variance, ts = E', '')
p38 = plot_heatmap_var(bEHLH_df_E_impact, c('bEH', 'LH'), 'Impact variance, ts = E', '')

p39 = plot_heatmap_var(bEHmuE_df_B_impact, c('bEH', 'muE'), 'Impact variance, ts = B', '')
p40 = plot_heatmap_var(bEHmuH_df_B_impact, c('bEH', 'muH'), 'Impact variance, ts = B', '')
p41 = plot_heatmap_var(bEHLH_df_B_impact, c('bEH', 'LH'), 'Impact variance, ts = B', '')

p42 = plot_heatmap_var(bEHmuE_df_H_impact, c('bEH', 'muE'), 'Impact variance, ts = H', '')
p43 = plot_heatmap_var(bEHmuH_df_H_impact, c('bEH', 'muH'), 'Impact variance, ts = H', '')
p44 = plot_heatmap_var(bEHLH_df_H_impact, c('bEH', 'LH'), 'Impact variance, ts = H', '')

p45 = plot_heatmap_var(bEHmuE_df_A_impact, c('bEH', 'muE'), 'Impact variance, ts = A', '')
p46 = plot_heatmap_var(bEHmuH_df_A_impact, c('bEH', 'muH'), 'Impact variance, ts = A', '')
p47 = plot_heatmap_var(bEHLH_df_A_impact, c('bEH', 'LH'), 'Impact variance, ts = A', '')

svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/viplot.svg', height=12, width = 20)
viplot = grid.arrange(p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,nrow=4,ncol=3)
dev.off()
"""

#And now recreate plots with varying ΛA
#1. change initial values of ΛA
p_E[:,2] .= rand(Uniform(0.000001, 1.), size(p_E)[1])
p_B[:,2] .= rand(Uniform(0.000001, 1.), size(p_E)[1])
p_H[:,2] .= rand(Uniform(0.000001, 1.), size(p_E)[1])
p_A[:,2] .= rand(Uniform(0.000001, 1.), size(p_E)[1])

#2. Re-run models
dat_E_2 = model_run(p_E, unboundeds)
dat_B_2 = model_run(p_B, unboundeds)
dat_H_2 = model_run(p_H, unboundeds)
dat_A_2 = model_run(p_A, unboundeds)

#3. get impact and summarise for parameter bins
impact_E_2 = [ifelse(!(dat_E_2[i]==0.), #can't have a 0 in denominator, and if 0 wouldn't expect any impact anyway
                   1 - dat_E_noLA[i]/dat_E_2[i], #will get %decrease and %increase here, will maybe cut out later
                   0) for i in 1:size(p_E)[1]]
impact_B_2 = [ifelse(!(dat_B_2[i]==0.),
                  1 - dat_B_noLA[i]/dat_B_2[i],
                  0) for i in 1:size(p_E)[1]]
impact_H_2 = [ifelse(!(dat_H_2[i]==0.),
                 1 - dat_H_noLA[i]/dat_H_2[i],
                 0) for i in 1:size(p_E)[1]]
impact_A_2 = [ifelse(!(dat_A_2[i]==0.),
                1 - dat_A_noLA[i]/dat_A_2[i],
                0) for i in 1:size(p_E)[1]]

@rput impact_B_2 impact_H_2 impact_A_2 impact_E_2 p_E p_H p_A p_B
R"""
impact_bEHLA_E = get_mean_var(lower_bin, p_E[, c(11, 2)], impact_E_2)
impact_bEHLA_H = get_mean_var(lower_bin, p_H[, c(11, 2)], impact_H_2)
impact_bEHLA_A = get_mean_var(lower_bin, p_A[, c(11, 2)], impact_A_2)
impact_bEHLA_B = get_mean_var(lower_bin, p_B[, c(11, 2)], impact_B_2)

p48 = plot_heatmap(impact_bEHLA_E, c('bEH', 'LA'), 'Mean impact, ts = E', '')
p49 = plot_heatmap(impact_bEHLA_B, c('bEH', 'LA'), 'Mean impact, ts = B', '')
p50 = plot_heatmap(impact_bEHLA_H, c('bEH', 'LA'), 'Mean impact, ts = H', '')
p51 = plot_heatmap(impact_bEHLA_A, c('bEH', 'LA'), 'Mean impact, ts = A', '')

p52 = plot_heatmap_var(impact_bEHLA_E, c('bEH', 'LA'), 'Impact variance, ts = E', '')
p53 = plot_heatmap_var(impact_bEHLA_B, c('bEH', 'LA'), 'Impact variance, ts = B', '')
p54 = plot_heatmap_var(impact_bEHLA_H, c('bEH', 'LA'), 'Impact variance, ts = H', '')
p55 = plot_heatmap_var(impact_bEHLA_A, c('bEH', 'LA'), 'Impact variance, ts = A', '')
"""

R"""
svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miLAplot.svg', height=12, width = 7)
miLAplot = grid.arrange(p48, p49, p50, p51, nrow=4)
dev.off()

svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/viLAplot.svg', height=12, width = 7)
viLAplot = grid.arrange(p52, p53, p54, p55, nrow=4)
dev.off()
"""
