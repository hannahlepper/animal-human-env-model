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

#plot each transmission scenario
u0 = [0.0; 0.0; 0.0]
tspan = (0.0, 1000.)
#       LH   LA   gH     gA     bHH         bAA         bHA         bAH         bAE         bEA         bEH         bHE
p_B =  [0.1, 0.1, 0.001, 0.001, 0.1,        0.1,        0.001,      0.1,        0.1,        0.01,       0.01,       0.1,        0.1, 0.1, 0.2]
p_Bd = [0.1, 0.1, 0.001, 0.001, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.07432092, 0.1, 0.1, 0.2]
p_E =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.001,      0.001,      0.001,      0.1420501,  0.1420501,  0.1420501,  0.1420501,  0.1, 0.1, 0.2]
p_A =  [0.1, 0.1, 0.001, 0.001, 0.001,      0.2019663,  0.001,      0.2019663,  0.2019663,  0.001,      0.001,      0.001,      0.1, 0.1, 0.2]
p_H =  [0.1, 0.1, 0.001, 0.001, 0.2019663,  0.001,      0.2019663,  0.001,      0.001,      0.001,      0.001,      0.2019663,  0.1, 0.1, 0.2]
p = [p_B, p_Bd, p_E, p_A, p_H]

R = [solve(ODEProblem(unboundeds, u0, tspan, p[i]))(200)[j] for i in 1:5, j in 1:3]

R"""
library(ggplot2)
library(readr)
RVals <- read_csv("M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/Manuscript\\/Fig\ 1\\/RVals.csv")
svg('M:/Project folders/Model env compartment/Plots/Manuscript/Fig 1/RVals.svg', height = 3, width = 5)
ggplot(RVals, aes(TS,Val,fill=R)) + geom_col(position = "dodge")
#dev.off()
"""
R"""dev.off()"""
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

#σ=1.
#pv = vp(ts(lnp(σ, 0.01), lnp(σ, 0.001),
#            lnp(σ, 0.14), lnp(σ, 0.001)), #βEH
#          lnp(σ, 0.1), #ΛH
#          lnp(σ, 0.2), #μE
#          lnp(σ, 0.1)) #μH

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
p_initial_Bd = copy(p_initial)
p_initial_H = copy(p_initial)
p_initial_A = copy(p_initial)
p_initial_E = copy(p_initial)
p_initial = [p_initial_B, p_initial_Bd, p_initial_H, p_initial_A, p_initial_E]

#Parameters that have sampling distributions
#p_initial_B[:,11] .= rand(LogNormal(log(pv.βEH.B.μ)+pv.βEH.B.σ, pv.βEH.B.σ), N)
#p_initial_H[:,11] .= rand(LogNormal(log(pv.βEH.H.μ)+pv.βEH.H.σ, pv.βEH.H.σ), N)
#p_initial_A[:,11] .= rand(LogNormal(log(pv.βEH.A.μ)+pv.βEH.A.σ, pv.βEH.A.σ), N)
#p_initial_E[:,11] .= rand(LogNormal(log(pv.βEH.E.μ)+pv.βEH.E.σ, pv.βEH.E.σ), N)
p_initial_B[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_Bd[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_H[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_A[:,11] .= rand(Uniform(0.000001, 1.), N)
p_initial_E[:,11] .= rand(Uniform(0.000001, 1.), N)

#Parameters that are fixed
#βHH, βAA, βHA, βAH, βAE, βEA, βHE
p_initial_B[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).B for x in [4,5,6,7,8,9,10]]'
p_initial_Bd[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).Bd for x in [4,5,6,7,8,9,10]]'
p_initial_H[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).H for x in [4,5,6,7,8,9,10]]'
p_initial_A[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).A for x in [4,5,6,7,8,9,10]]'
p_initial_E[:,[5,6,7,8,9,10,12]] .= [getfield(pf,x).E for x in [4,5,6,7,8,9,10]]'

#get rid of negative numbers, or values over 1.5
function keep_ps(p)
    sets_LH = p[:,1] .< 1.
    sets_muH = 0.002 .< p[:,13] .< 1.
    sets_muE = p[:,15] .<1.5
    keep = findall(sets_LH .& sets_muH .& sets_muE)
    return p[keep,:]
end
p_B = keep_ps(p_initial_B)
p_Bd = keep_ps(p_initial_Bd)
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
# 7. = fixed bEH, fixed LA = 0.1

#2. varying bEH, fixed LA = 0.0
p_E_2 = copy(p_E)
p_B_2 = copy(p_B)
p_Bd_2 = copy(p_Bd)
p_H_2 = copy(p_H)
p_A_2 = copy(p_A)
p_B_2[:,2] .= 0
p_Bd_2[:,2] .= 0
p_H_2[:,2] .= 0
p_A_2[:,2] .= 0
p_E_2[:,2] .= 0


#3. varying LA, varying bEH
p_E_3 = copy(p_E)
p_B_3 = copy(p_B)
p_Bd_3 = copy(p_Bd)
p_H_3 = copy(p_H)
p_A_3 = copy(p_A)
p_E_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_E)[1])
p_B_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_B)[1])
p_Bd_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_Bd)[1])
p_H_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_H)[1])
p_A_3[:,2] .= rand(Uniform(0.000001, 1.), size(p_A)[1])

#4. fixed bEH, varying LA
p_E_4 = copy(p_E_3)
p_B_4 = copy(p_B_3)
p_Bd_4 = copy(p_Bd_3)
p_H_4 = copy(p_H_3)
p_A_4 = copy(p_A_3)
p_E_4[:,11] .= 0.1420501
p_B_4[:,11] .= 0.01
p_Bd_4[:,11] .= 0.01
p_H_4[:,11] .= 0.001
p_A_4[:,11] .= 0.001

#5. fixed bEH = 0, fixed LA = 0.1
p_E_5 = copy(p_E)
p_B_5 = copy(p_B)
p_Bd_5 = copy(p_Bd)
p_H_5 = copy(p_H)
p_A_5 = copy(p_A)
p_E_5[:,11] .= 0
p_B_5[:,11] .= 0
p_Bd_5[:,11] .= 0
p_H_5[:,11] .= 0
p_A_5[:,11] .= 0

#6. fixed bEH, fixed LA = 0
p_E_6 = copy(p_E_4)
p_B_6 = copy(p_B_4)
p_Bd_6 = copy(p_Bd_4)
p_H_6 = copy(p_H_4)
p_A_6 = copy(p_A_4)
p_B_6[:,2] .= 0
p_Bd_6[:,2] .= 0
p_H_6[:,2] .= 0
p_A_6[:,2] .= 0
p_E_6[:,2] .= 0

#7. fixed bEH = ts, fixed LA = 0.1
p_E_7 = copy(p_E_6)
p_B_7 = copy(p_B_6)
p_Bd_7 = copy(p_Bd_6)
p_H_7 = copy(p_H_6)
p_A_7 = copy(p_A_6)
p_Bd_7[:,2] .= 0.1
p_H_7[:,2] .= 0.1
p_A_7[:,2] .= 0.1
p_E_7[:,2] .= 0.1
