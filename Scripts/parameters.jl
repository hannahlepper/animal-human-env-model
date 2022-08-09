
using Distributions
#using Plots
#using RCall

include("optim_parameters_HAE.jl") 

#initial parameter values
#       LH   LA   gH     gA     bHH         bAA         bHA         bAH         bAE         bEA         bEH         bHE         muH  muA  muE
p_Bd = [0.1, 0.1, 0.001, 0.001, b_bounded_bal, b_bounded_bal, b_bounded_bal, b_bounded_bal, b_bounded_bal, b_bounded_bal, b_bounded_bal, b_bounded_bal, 0.2, 0.4, 0.29]
p_E =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.001,      0.001,      0.001,      b_bounded_env,  b_bounded_env,  b_bounded_env,  b_bounded_env,  0.2, 0.4, 0.29]
p_A =  [0.1, 0.1, 0.001, 0.001, 0.001,      b_bounded_anim,  0.001,      b_bounded_anim,  b_bounded_anim,  0.001,      0.001,      0.001,      0.2, 0.4, 0.29]
p_H =  [0.1, 0.1, 0.001, 0.001, b_bounded_hum,  0.001,      b_bounded_hum,  0.001,      0.001,      0.001,      0.001,      b_bounded_hum,  0.2, 0.4, 0.29]
p = [p_Bd, p_E, p_A, p_H]

# R"""
# library(ggplot2)
# library(readr)
# RVals <- read_csv("M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/Manuscript\\/Fig\ 1\\/RVals.csv")
# svg('M:/Project folders/Model env compartment/Plots/Manuscript/Fig 1/RVals.svg', height = 3, width = 5)
# ggplot(RVals, aes(TS,Val,fill=R)) + geom_col(position = "dodge")
# #dev.off()
# """
# R"""dev.off()"""

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
    Bd
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

#beta distribution with mean = parameter selected in transmission scenarios, and var = mean/20
#ΛH 
LHmu = 0.1 
LHvar = LHmu/20
LH_alpha = LHmu * ((LHmu * (1 - LHmu))/LHvar - 1)
LH_beta = (1 - LHmu) * ((LHmu * (1 - LHmu))/LHvar - 1)

#μE 
mEmu = 0.29 
mEvar = mEmu/20
mE_alpha = mEmu * ((mEmu * (1 - mEmu))/mEvar - 1)
mE_beta = (1 - mEmu) * ((mEmu * (1 - mEmu))/mEvar - 1)

#μH 
mHmu = 0.118 
mHvar = mHmu/20
mH_alpha = mHmu * ((mHmu * (1 - mHmu))/mHvar - 1)
mH_beta = (1 - mHmu) * ((mHmu * (1 - mHmu))/mHvar - 1)

# order reminder - p_Bd, p_H, p_E, p_A
bEH_mu = [b_unbounded_bal, 0.001, b_unbounded_env, 0.001]
bEH_var = bEH_mu./20

bEH_alpha = map(i -> bEH_mu[i] * ((bEH_mu[i] * (1 - bEH_mu[i]))/bEH_var[i] - 1), 1:4)
bEH_beta = map(i -> (1 - bEH_mu[i]) * ((bEH_mu[i] * (1 - bEH_mu[i]))/bEH_var[i] - 1), 1:4)

pv = vp(ts(bp(bEH_alpha[1],bEH_beta[1]), 
           bp(bEH_alpha[2],bEH_beta[2]), 
           bp(bEH_alpha[3], bEH_beta[3]),
           bp(bEH_alpha[4], bEH_beta[4])), #βEH
        bp(LH_alpha, LH_beta), #ΛH
        bp(mE_alpha, mE_beta), #μE
        bp(mH_alpha, mH_beta)) #μH
bEH_ts_1 = bEH_mu


bEH_mu = [b_bounded_bal, 0.001, b_bounded_env, 0.001]
bEH_var = bEH_mu./20

bEH_alpha = map(i -> bEH_mu[i] * ((bEH_mu[i] * (1 - bEH_mu[i]))/bEH_var[i] - 1), 1:4)
bEH_beta = map(i -> (1 - bEH_mu[i]) * ((bEH_mu[i] * (1 - bEH_mu[i]))/bEH_var[i] - 1), 1:4)

pv2 = vp(ts(bp(bEH_alpha[1],bEH_beta[1]), 
           bp(bEH_alpha[2],bEH_beta[2]), 
           bp(bEH_alpha[3], bEH_beta[3]),
           bp(bEH_alpha[4], bEH_beta[4])), #βEH
        bp(LH_alpha, LH_beta), #ΛH
        bp(mE_alpha, mE_beta), #μE
        bp(mH_alpha, mH_beta)) #μH
bEH_ts_2 = bEH_mu

bEH_ts_3 = [b_orig_bal, 0.001, b_orig_env, 0.001]

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

#Unbounded
pf = f_p(0.1,0.001,0.001, #1 = ΛA, 2 = γH, 3 = γA

         ts(b_unbounded_bal, b_unbounded_hum,           0.001,            0.001), #4 = βHH
         ts(b_unbounded_bal,           0.001,           0.001, b_unbounded_anim), #5 = βAA
         ts(b_unbounded_bal, b_unbounded_hum,           0.001,            0.001), #6 = βHA
         ts(b_unbounded_bal,           0.001,           0.001, b_unbounded_anim), #βAH
         ts(b_unbounded_bal,           0.001, b_unbounded_env, b_unbounded_anim), #βAE
         ts(b_unbounded_bal,           0.001, b_unbounded_env,            0.001), #βEA
         ts(b_unbounded_bal, b_unbounded_hum, b_unbounded_env,            0.001), #βHE
         mHmu*2) #μA

#Bounded
pf2 = f_p(0.1,0.001,0.001, #ΛA, γH, γA
         ts(b_bounded_bal, b_bounded_hum,         0.001,          0.001), #4 = βHH
         ts(b_bounded_bal,         0.001,         0.001, b_bounded_anim), #5 = βAA
         ts(b_bounded_bal, b_bounded_hum,         0.001,          0.001), #6 = βHA
         ts(b_bounded_bal,         0.001,         0.001, b_bounded_anim), #βAH
         ts(b_bounded_bal,         0.001, b_bounded_env, b_bounded_anim), #βAE
         ts(b_bounded_bal,         0.001, b_bounded_env,          0.001), #βEA
         ts(b_bounded_bal, b_bounded_hum, b_bounded_env,          0.001), #βHE
         mHmu*2) #μA

#Original
pf3 = f_p(0.1,0.001,0.001, #ΛA, γH, γA
         ts(b_orig_bal, b_orig_hum,       0.001,       0.001), #4 = βHH
         ts(b_orig_bal,       0.001,      0.001, b_orig_anim), #5 = βAA
         ts(b_orig_bal, b_orig_hum,       0.001,       0.001), #6 = βHA
         ts(b_orig_bal,       0.001,      0.001, b_orig_anim), #βAH
         ts(b_orig_bal,       0.001, b_orig_env, b_orig_anim), #βAE
         ts(b_orig_bal,       0.001, b_orig_env,       0.001), #βEA
         ts(b_orig_bal,  b_orig_hum, b_orig_env,       0.001), #βHE
         mHmu*2) #μA

#set up parameter set that never changes first
N = 2000000 
p_uncertainty = zeros(N, 3)
p_uncertainty2 = zeros(N, 3)
p_uncertainty3 = zeros(N, 3)

#Parameters varying to account for uncertainty
using Random
Random.seed!(123)

p_uncertainty[:,1] .= rand(Beta(pv.ΛH.α, pv.ΛH.β), N)
p_uncertainty[:,2] .= rand(Beta(pv.μE.α, pv.μE.β), N)
p_uncertainty[:,3] .= rand(Beta(pv.μH.α, pv.μH.β), N)

p_uncertainty2[:,1] .= rand(Beta(pv2.ΛH.α, pv2.ΛH.β), N)
p_uncertainty2[:,2] .= rand(Beta(pv2.μE.α, pv2.μE.β), N)
p_uncertainty2[:,3] .= rand(Beta(pv2.μH.α, pv2.μH.β), N)

p_uncertainty3[:,1] .= rand(Beta(pv2.ΛH.α, pv2.ΛH.β), N)
p_uncertainty3[:,2] .= rand(Beta(pv2.μE.α, pv2.μE.β), N)
p_uncertainty3[:,3] .= rand(Beta(pv2.μH.α, pv2.μH.β), N)


function col_edit(mat, col_index, replacement)
    mat[:,col_index] .= replacement
    return mat
end

#Experiment arrays

bEH_unif = rand(Uniform(0.000001, 1.), N)
bEH_experiments_unbounded = [
    bEH_ts_1,  bEH_ts_1,   bEH_ts_1,
    0.1,    0.5,  0,
    bEH_unif,    bEH_unif,    bEH_unif,   bEH_unif, 0, 0, bEH_unif
]

bEH_experiments_bounded = [
    bEH_ts_2,   bEH_ts_2,   bEH_ts_2,
    0.1,    0.5,   0,
    bEH_unif,    bEH_unif,    bEH_unif,    bEH_unif, 0, 0, bEH_unif
] 

bEH_experiments_orig = [0,0,0,0,0,0,0,0,0,0,0,0,0]

LA_unif = rand(Uniform(0.000001, 1.), N)
LA_experiments = [
    0.1, 0.0, 0.5, 0.1, 0.1, 0.1, 0., LA_unif, 0., LA_unif, 0, 0, 0.1
]

#Index for runs
# 1. = fixed bEH, LA fixed to 0.1.
# 2. = fixed bEH, LA fixed to 0.0.
# 3. = fixed bEH, LA fixed to 0.5. mistake in previous code here - p_2 has varying bEH and varying lA?
# 4. = bEH fixed to 0.1, fixed LA = 0.1. mistake in previous code - fixed bEH, varying LA
# 5. = bEH fixed to 0.5, fixed LA = 0.1. mistake in previous code - fixed bEH = 0, fixed LA = 0.1
# 6. = bEH fixed to 0, fixed LA = 0.1 - also needed for  fig 3.C. mistake in previous code - fixed bEH, fixed LA = 0
# 7. = varying bEH, low bHA, LA fixed to 0. mistake in previous code - fixed bEH = ts, fixed LA = 0.1
# 8. = varying bEH, low bHA, varying LA. previous code - varying bEH, low bHA, LA fixed to 0.
# 9. = varying bEH, high bHA, LA fixed to 0. previous code - varying bEH, low bHA, varying LA.
# 10. = varying bEH, high bHA, varying LA. previous code - bEH and LA set to 0
# 11. = bEH fixed to 0, LA fixed to 0, high/normal bHA
# 12. = bEH fixed to 0, LA fixed to 0, low bHA
# 13. = bEH varying, LA fixed to 1, normal bHA

#Figure 3C data needed. Make separately so that
# original model - all environmnetal parameters set to 0, otherwise following baseline parameter values
# original model, with LA set to 0


#Function for getting parameters for a model run

function get_params(transmission_scenario, experiment_num, n_sets,
    fixed_pars, uncertainty_pars, 
    bEH_exp, LA_exp)
    
    #Order of parameters in model function: 
    #ΛH, ΛA, γH, γA, βHH, βAA, βHA, βAH, βAE, βEA, βEH, βHE, μH, μA, μE
    p_mat = zeros(n_sets, 15)

    #Parameters that are the same in every run: ΛA, γH, γA, μA
    p_mat[:,[2,3,4,14]] .= [fixed_pars.ΛA, fixed_pars.γH, fixed_pars.γA, fixed_pars.μA]'

    #Parameters that do not vary in this experiment
    #1. Transmission scenario related parameters
    #βHH, βAA, βHA, βAH, βAE, βEA, βHE
    p_mat[:, [5,6,7,8,9,10,12]] .= [getfield(fixed_pars.βHH, transmission_scenario), 
                                    getfield(fixed_pars.βAA, transmission_scenario), 
                                    getfield(fixed_pars.βHA, transmission_scenario), 
                                    getfield(fixed_pars.βAH, transmission_scenario),
                                    getfield(fixed_pars.βAE, transmission_scenario), 
                                    getfield(fixed_pars.βEA, transmission_scenario), 
                                    getfield(fixed_pars.βHE, transmission_scenario)]'

    #2. Experiment-related fixed parameters
    #bEH
    if length(bEH_exp[experiment_num]) == 4
        p_mat[:,11] .= bEH_exp[experiment_num][transmission_scenario]
    elseif length(bEH_exp[experiment_num]) > 4
        p_mat[:,11] = bEH_exp[experiment_num][1:n_sets]
    elseif length(bEH_exp[experiment_num]) == 1
        p_mat[:,11] .= bEH_exp[experiment_num]
    end

    #LA
    if length(LA_exp[experiment_num]) > 1
        p_mat[:,2] = LA_exp[experiment_num][1:n_sets]
    else
        p_mat[:,2] .= LA_exp[experiment_num]
    end 

    #bHA
    if in([7,8, 12]).(experiment_num)
        p_mat[:,7] .= p_mat[:,7]./100
    end 

    #Uncertainty parameters
    #ΛH, μE, μH
    p_mat[:, [1,15,13]] .= uncertainty_pars[1:n_sets, :]

    return p_mat

end

get_params(1, 13, 2, pf3, p_uncertainty3, bEH_experiments_orig, LA_experiments)

#11. Original model

# 1  2  3  4  5   6   7   8   9   10  11  12  13  14  15
# LH LA gH gA bHH bAA bHA bAH bAE bEA bEH bHE muH muA muE
     
#p_orig = get_params(1, 1, 2000000, pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments)
#p_orig = col_edit(p_orig, [3,4,9,10,11,12,15], 0)

#p_orig_int = deepcopy(p_orig)
#p_orig_int = col_edit(p_orig_int, 2, 0)

# using JLD2
# @save "/mnt/d/results_workspace_p_orig.jld2" p_orig
# @save "/mnt/d/results_workspace_p_orig_int.jld2" p_orig_int
