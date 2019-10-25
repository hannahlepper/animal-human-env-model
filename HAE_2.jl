using DifferentialEquations
using Distributions
using Plots
using CSV
using Tables
using ORCA #needed for plotlyjs

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
    tspan = (0.,2000.)
    for i in re_run
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

#get indexes of p rows where this is true and calculate %
function get_perc_target(bins, binsize, p, dat, p1, p2)
    dims_out = size(bins)[1]
    dims_in = size(dat)[1]
    perc_target = zeros(dims_out, dims_out)
    for i in 1:dims_out
        for j in 1:dims_out
            sum_success=0
            n_sets=0
            for k in 1:dims_in
                if p[k,p1] > lower_bin[i] && p[k,p1] <  lower_bin[i] + binsize #bEH
                    if p[k,p2] > lower_bin[j] && p[k,p2] < lower_bin[j] + binsize #muE
                        n_sets += 1
                        sum_success += dat[k]
                    end
                end
            end
            perc_target[j,i] = ifelse(sum_success > 0 && n_sets > 0, sum_success/n_sets, 0)
        end
    end
    return perc_target
end
#run on small thing to get function ready
get_perc_target(lower_bin, 0.05, p_E[1:100,:], n_target[1:100],1,2)


# Conclusion 1: realistic RHs are attainable for environmental transmission scenarios
@time βEHμE_mat_E = get_perc_target(lower_bin,0.05, p_E, n_target, 11, 15)
@time βEHμH_mat_E = get_perc_target(lower_bin,0.05, p_E, n_target, 11, 13)
@time βEHΛH_mat_E = get_perc_target(lower_bin,0.05, p_E, n_target, 11, 1)

@time βEHμE_mat_B = get_perc_target(lower_bin,0.05, p_B, n_target, 11, 15)
@time βEHμH_mat_B = get_perc_target(lower_bin,0.05, p_B, n_target, 11, 13)
@time βEHΛH_mat_B = get_perc_target(lower_bin,0.05, p_B, n_target, 11, 1)

@time βEHμE_mat_H = get_perc_target(lower_bin,0.05, p_H, n_target, 11, 15)
@time βEHμH_mat_H = get_perc_target(lower_bin,0.05, p_H, n_target, 11, 13)
@time βEHΛH_mat_H = get_perc_target(lower_bin,0.05, p_H, n_target, 11, 1)

@time βEHμE_mat_A = get_perc_target(lower_bin,0.05, p_A, n_target, 11, 15)
@time βEHμH_mat_A = get_perc_target(lower_bin,0.05, p_A, n_target, 11, 13)
@time βEHΛH_mat_A = get_perc_target(lower_bin,0.05, p_A, n_target, 11, 1)


#plotting!
#Env
function ph(mat, x, y)
    heatmap(lower_bin, lower_bin, mat, fillcolor = :fire,
            xlabel = x, ylabel = y)
end

p1 = ph(βEHμE_mat_E, "βEH", "μE")
p2 = ph(βEHΛH_mat_E, "βEH", "ΛH")
p3 = ph(βEHμH_mat_E, "βEH", "μH")

#Base
p4 = ph(βEHμE_mat_B, "βEH", "μE")
p5 = ph(βEHΛH_mat_B, "βEH", "ΛH")
p6 = ph(βEHμH_mat_B, "βEH", "μH")

#Human
p7 = ph(βEHμE_mat_H, "βEH", "μE")
p8 = ph(βEHΛH_mat_H, "βEH", "ΛH")
p9 = ph(βEHμH_mat_H, "βEH", "μH")

#Animal
p10 = ph(βEHμE_mat_A, "βEH","μE")
p11 = ph(βEHμH_mat_A, "βEH","ΛH")
p12 = ph(βEHΛH_mat_A, "βEH","μH")

p_bEH_target = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,layout = (4,3))

#none of this is working but I want to move on
PlotlyJS.savefig(p1.o, "M:/Project folders/Model env compartment/Plots/EbEHmuEheat.svg")
PlotlyJS.savefig(p2.o, "M:/Project folders/Model env compartment/Plots/EbEHmuHheat.svg")
PlotlyJS.savefig(p3.o, "M:/Project folders/Model env compartment/Plots/EbEHLHheat.svg")
plotlyjs()
pgfplots()
savefig(p_bEH_target, "M:/Project folders/Model env compartment/Plots/bEHTargetheat.svg")
p_bEH_target.o

#Put into CSVs so I don't have to run again unless I want to
CSV.write("M:/Project folders/Model env compartment/bEHmE.csv", Tables.table(βEHμE_mat))
CSV.write("M:/Project folders/Model env compartment/bEHmH.csv", Tables.table(βEHμH_mat))
CSV.write("M:/Project folders/Model env compartment/bEHLH.csv", Tables.table(βEHΛH_mat))


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

@time βEHμE_low_impact_B = get_perc_target(lower_bin,0.05, p_B, n_lowimpact_B, 11, 15)
@time βEHμH_low_impact_B = get_perc_target(lower_bin,0.05, p_B, n_lowimpact_B, 11, 13)
@time βEHΛH_low_impact_B = get_perc_target(lower_bin,0.05, p_B, n_lowimpact_B, 11, 1)

@time βEHμE_low_impact_H = get_perc_target(lower_bin,0.05, p_H, n_lowimpact_H, 11, 15)
@time βEHμH_low_impact_H = get_perc_target(lower_bin,0.05, p_H, n_lowimpact_H, 11, 13)
@time βEHΛH_low_impact_H = get_perc_target(lower_bin,0.05, p_H, n_lowimpact_H, 11, 1)

@time βEHμE_low_impact_A = get_perc_target(lower_bin,0.05, p_A, n_lowimpact_A, 11, 15)
@time βEHμH_low_impact_A = get_perc_target(lower_bin,0.05, p_A, n_lowimpact_A, 11, 13)
@time βEHΛH_low_impact_A = get_perc_target(lower_bin,0.05, p_A, n_lowimpact_A, 11, 1)

@time βEHμE_low_impact_E = get_perc_target(lower_bin,0.05, p_E, n_lowimpact_E, 11, 15)
@time βEHμH_low_impact_E = get_perc_target(lower_bin,0.05, p_E, n_lowimpact_E, 11, 13)
@time βEHΛH_low_impact_E = get_perc_target(lower_bin,0.05, p_E, n_lowimpact_E, 11, 1)

p13 = ph(βEHμE_low_impact_B, "βEH","μE")
p14 = ph(βEHμH_low_impact_B, "βEH","μH")
p15 = ph(βEHΛH_low_impact_B, "βEH","ΛH")
p16 = ph(βEHμE_low_impact_H, "βEH","μE")
p17 = ph(βEHμH_low_impact_H, "βEH","μH")
p18 = ph(βEHΛH_low_impact_H, "βEH","ΛH")
p19 = ph(βEHμE_low_impact_A, "βEH","μE")
p20 = ph(βEHμH_low_impact_A, "βEH","μH")
p21 = ph(βEHΛH_low_impact_A, "βEH","ΛH")
p22 = ph(βEHμE_low_impact_E, "βEH","μE")
p23 = ph(βEHμH_low_impact_E, "βEH","μH")
p24 = ph(βEHΛH_low_impact_E, "βEH","ΛH")

βEH_low_impact = plot(p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, layout = (4,3))



PlotlyJS.savefig(p4.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuEheat2.svg")

PlotlyJS.savefig(p5.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuHheat2.svg")

PlotlyJS.savefig(p6.o,
        "M:/Project folders/Model env compartment/Plots/bEHLHheat2.svg")

#Final conclusion - increasing βEH in many cases reduces the impacts of increased LH
# need function for averaging the bins

function bin_av(dat, lower_bin, p, rows, cols)
    #1. make array indicating the position of parameters in rows and cols of final matrix.
    #Bear in mind, in this function, all values higher than 1.0 go in the final bin.
    row_bin_index = [findlast(p[i,rows] .> lower_bin) for i in 1:size(p)[1]]
    col_bin_index = [findlast(p[i,cols] .> lower_bin) for i in 1:size(p)[1]]
    #2. make the matrix
    mat = [mean(dat[findall((row_bin_index.==r) .& (col_bin_index.==c))]) for r in 1:size(lower_bin)[1], c in 1:size(lower_bin)[1]]
    return mat
end

function bin_var(dat, lower_bin, p, rows, cols)
    #1. make array indicating the position of parameters in rows and cols of final matrix.
    #Bear in mind, in this function, all values higher than 1.0 go in the final bin.
    row_bin_index = [findlast(p[i,rows] .> lower_bin) for i in 1:size(p)[1]]
    col_bin_index = [findlast(p[i,cols] .> lower_bin) for i in 1:size(p)[1]]
    #2. make the matrix
    mat = [var(dat[findall((row_bin_index.==r) .& (col_bin_index.==c))]) for r in 1:size(lower_bin)[1], c in 1:size(lower_bin)[1]]
    return mat
end

impact_βEHΛH_E = bin_av(impact_E, lower_bin, p_E, 1, 11)
impact_βEHμH_E = bin_av(impact_E, lower_bin, p_E, 1, 13)
impact_βEHμE_E = bin_av(impact_E, lower_bin, p_E, 1, 15)

impact_βEHΛH_H = bin_av(impact_H, lower_bin, p_H, 1, 11)
impact_βEHμH_H = bin_av(impact_H, lower_bin, p_H, 1, 13)
impact_βEHμE_H = bin_av(impact_H, lower_bin, p_H, 1, 15)

impact_βEHΛH_A = bin_av(impact_A, lower_bin, p_A, 1, 11)
impact_βEHμH_A = bin_av(impact_A, lower_bin, p_A, 1, 13)
impact_βEHμE_A = bin_av(impact_A, lower_bin, p_A, 1, 15)

impact_βEHΛH_B = bin_av(impact_B, lower_bin, p_B, 1, 11)
impact_βEHμH_B = bin_av(impact_B, lower_bin, p_B, 1, 13)
impact_βEHμE_B = bin_av(impact_B, lower_bin, p_B, 1, 15)

p25 = ph(impact_βEHμE_B, "βEH","μE")
p26 = ph(impact_βEHμH_B, "βEH","μH")
p27 = ph(impact_βEHΛH_B, "βEH","ΛH")
p28 = ph(impact_βEHμE_H, "βEH","μE")
p29 = ph(impact_βEHμH_H, "βEH","μH")
p30 = ph(impact_βEHΛH_H, "βEH","ΛH")
p31 = ph(impact_βEHμE_A, "βEH","μE")
p32 = ph(impact_βEHμH_A, "βEH","μH")
p33 = ph(impact_βEHΛH_A, "βEH","ΛH")
p34 = ph(impact_βEHμE_E, "βEH","μE")
p35 = ph(impact_βEHμH_E, "βEH","μH")
p36 = ph(impact_βEHΛH_E, "βEH","ΛH")
