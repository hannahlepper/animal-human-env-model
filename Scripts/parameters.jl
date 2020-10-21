
using Distributions
#using Plots
#using RCall

#initial parameter values
#       LH   LA   gH     gA     bHH         bAA         bHA         bAH         bAE         bEA         bEH         bHE         muH  muA  muE
p_B =  [0.1, 0.1, 0.001, 0.001, 0.1,        0.1,        0.001,      0.1,        0.1,        0.01,       0.01,       0.1,        0.1, 0.1, 0.2]
p_Bd = [0.1, 0.1, 0.001, 0.001, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.1, 0.1, 0.2]
p_E =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.001,      0.001,      0.001,      0.1420501,  0.1420501,  0.1420501,  0.1420501,  0.1, 0.1, 0.2]
p_A =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.2019663,  0.001,      0.2019663,  0.2019663,  0.001,      0.001,      0.001,      0.1, 0.1, 0.2]
p_H =  [0.1, 0.1, 0.001, 0.001, 0.2019663,  0.001,      0.2019663,  0.001,      0.001,      0.001,      0.001,      0.2019663,  0.1, 0.1, 0.2]
p = [p_B, p_Bd, p_E, p_A, p_H]

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
    B
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
#parameters estimated from R code
pv = vp(ts(bp(0.188, 18.612), bp(1.301625,16.21196), bp(0.01898, 18.96102),
          bp(2.295387, 13.86361), bp(0.01898, 18.96102)), #βEH
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
         ts(0.1,   0.07432092, 0.2019663,     0.001,     0.001), #βHH
         ts(0.1,   0.07432092,     0.001,     0.001, 0.2019663), #βAA
         ts(0.001, 0.07432092, 0.2019663,     0.001,     0.001), #βHA
         ts(0.1,   0.07432092,     0.001,     0.001, 0.2019663), #βAH
         ts(0.1,   0.07432092,     0.001, 0.1420501, 0.2019663), #βAE
         ts(0.01,  0.07432092,     0.001, 0.1420501,     0.001), #βEA
         ts(0.1,   0.07432092, 0.2019663, 0.1420501,     0.001), #βHE
         0.1) #μA

#set up parameter set that never changes first
N = 2000000 
p_initial = zeros(N, 15)

using Random
Random.seed!(123)
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

p_initial
#Make parameter sets for the different transmission scenarios
p_1 = collect((map(x -> deepcopy(p_initial), 1:5)))

function col_edit(mat, col_index, replacement)
    mat[:,col_index] .= replacement
    return mat
end

#Parameters that have sampling distributions
p_1 = map(x -> col_edit(p_1[x], 11, rand(Uniform(0.000001, 1.), N)), 1:5)

#Parameters that are fixed
#βHH, βAA, βHA, βAH, βAE, βEA, βHE
#Order of transmission scenarios: B, Bd, H, E, A
p_1 = map(x -> col_edit(p_1[x], [5,6,7,8,9,10,12], [getfield(getfield(pf,p),x) for p in [4,5,6,7,8,9,10]]'), 1:5)

#get rid of negative numbers, or values over 1.5
function keep_ps(p)
    sets_LH = p[:,1] .< 1.
    sets_muH = 0.002 .< p[:,13] .< 1.
    sets_muE = p[:,15] .<1.5
    keep = findall(sets_LH .& sets_muH .& sets_muE)
    return p[keep,:]
end
p_1= map(x -> keep_ps(p_1[x]), 1:5)

#Index for runs
# Figure 2 data needed
# 1. = fixed bEH, LA fixed to 0.1.
# 2. = fixed bEH, LA fixed to 0.0.
# 3. = fixed bEH, LA fixed to 0.5.
# 4. = bEH fixed to 0.1, fixed LA = 0.1
# 5. = bEH fixed to 0.5, fixed LA = 0.1
# 6. = bEH fixed to 0, fixed LA = 0.1 - also needed for  fig 3.C

#Figure 3 A and B data needed
# 7. = varying bEH, low bHA, LA fixed to 0.
# 8. = varying bEH, low bHA, varying LA.
# 9. = varying bEH, high bHA, LA fixed to 0.
# 10. = varying bEH, high bHA, varying LA.

#Figure 3C data needed. Make separately so that
# 11. original model - all environmnetal parameters set to 0, otherwise following baseline parameter values
# 12. original model, with LA set to 0
# 13. varying bEH 0 -> 1, LA fixed to 0.1
# 14. varying bEH 0 -> 1, LA fixed to 0.0

#2. varying bEH, fixed LA = 0.0
p_2 = deepcopy(p_1)
p_2 = map(x -> col_edit(p_2[x], 2, 0), 1:5)

#3. varying LA, varying bEH
p_3 = deepcopy(p_2)
p_3 = map(x -> col_edit(p_3[x],2,rand(Uniform(0.000001, 1.), N)[1:size(p_3[1])[1]]), 1:5)

#4. fixed bEH, varying LA
p_4 = deepcopy(p_3)
bEHts = [0.1420501, 0.01, 0.01, 0.001, 0.001]
p_4 = map(x -> col_edit(p_4[x], 11, bEHts[x]), 1:5)

#5. fixed bEH = 0, fixed LA = 0.1
p_5 = deepcopy(p_1)
p_5 = map(x -> col_edit(p_5[x],11, 0), 1:5)

#6. fixed bEH, fixed LA = 0
p_6 = deepcopy(p_4)
p_6 = map(x -> col_edit(p_6[x],2,0), 1:5)

#7. fixed bEH = ts, fixed LA = 0.1
p_7 = deepcopy(p_6)
p_7 = map(x -> col_edit(p_7[x],2,0.1), 1:5)

#8.varying bEH, low bHA, LA fixed to 0.
p_8 = deepcopy(p_2)
p_8 = map(x -> col_edit(p_8[x],7, p_8[x][:,7]/100), 1:5)

#9. varying bEH, low bHA, varying LA.
p_9 = deepcopy(p_8)
p_9 = map(x -> col_edit(p_9[x],2,rand(Uniform(0.000001, 1.), N)[1:size(p_9[1])[1]]), 1:5)

#10. bEH and LA set to 0
p_10 = deepcopy(p_6)
p_10 = map(x -> col_edit(p_10[x], 11, 0), 1:5)

P = hcat(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
Pv = vec(map(x -> vec(x), [p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10]))

#10. Original model
p_orig = deepcopy(p_1[1])
p_orig = col_edit(p_orig, [3,4,9,10,11,12,15], 0)

#11. Original model - post intervention
p_orig_int = deepcopy(p_orig)
p_orig_int = col_edit(p_orig_int, 2, 0)

p_origs_v = vec(map(x -> vec(x), [p_orig, p_orig_int]))
