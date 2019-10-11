using DifferentialEquations
using Distributions
using Plots
using CSV
using Tables
using StatsBase

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
plotlyjs()
plot(sol, ylims=(0.,1.), yticks=0.:.1:1.)

#sample paramter space
N = 2000000
σ = 1.
#parameters of interest are bEH, LH, muE, muH.
p_initial = zeros(N, 15)

# #use log-normal for the environmental parameters - should explore possibility they are high
# MPβEH = 0.14
# p_initial[:,11].= rand(LogNormal(log(MPβEH),σ), N)
# MPμE = 0.2
# p_initial[:,15] .= rand(LogNormal(log(MPμE), σ), N)
# #use beta dist with mean 0.1 and var 0.1 for the human parameters, which I think are less likely to be >1
# μH_mean = 0.1
# min_var = μH_mean * (1 - μH_mean)
# μH_var = min_var - 0.001
# #needs to be less than min above for the equations below to work. Assumes one mode, no antimode.
# α_beta = ((1 - μH_mean)/μH_var^2 - 1/μH_mean) * μH_mean^2
# β_beta = α_beta * (1/μH_mean - 1)
# p_initial[:,13] .= rand(Beta(α_beta, β_beta), N)
# p_initial[:,1] .= rand(Beta(α_beta, β_beta), N)

#Alternative
σ = 1.
MPΛH = 0.1
p_initial[:,1] .= rand(LogNormal(log(MPΛH)+σ, σ), N)
MPβEH = 0.14
p_initial[:,11].= rand(LogNormal(log(MPβEH)+σ,σ), N)
MPμE = 0.2
p_initial[:,15] .= rand(LogNormal(log(MPμE)+σ, σ), N)
MPμH = 0.1
p_initial[:,13] .= rand(LogNormal(log(MPμH)+σ, σ), N)

histogram(p_initial[:,[1,11,13,15]],
          xticks = range(0, 1.5; step =0.1),
          xlims = (0.0, 1.5),
          label = ["ΛH" "βEH" "μH" "μE"])

#Rest of the parameters, on tha basis of the
#ΛA and μA
p_initial[:,[2,14]] .= 0.1
#γH, γA, βHH, βAA, βHA, βAH
p_initial[:,[3,4,5,6,7,8]] .= 0.001
#βEA, βAE, βHE
p_initial[:,[9,10,12]] .= 0.14

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
p = keep_ps(p_initial)

histogram(p[:,[1,11,13,15]],
          xticks = range(0, 1.5; step =0.1),
          xlims = (0.0, 1.5),
          label = ["ΛH" "βEH" "μH" "μE"])

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
model_run(p[:,1:5], unboundeds) #short run to get the function going
dat = model_run(p, unboundeds)

#Find number of runs where 0.65 < RH < 0.75
n_target = zeros(size(p)[1])
for i in 1:size(p)[1]
    if 0.65 < dat[i,2] < 0.75
        n_target[i] = 1
    else
        n_target[i] = 0
    end
end

#Bins for parameters
lower_bin = [0.:0.05:1;]
bin_N = size(lower_bin)[1]

#get βEH and μE combinations
βEHμE = zeros(bin_N^2,2)
βEHμE[:,1] = repeat(lower_bin; outer=bin_N)
βEHμE[:,2] = repeat(lower_bin; inner=bin_N)

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
#run on small thing first
get_perc_target(lower_bin, 0.05, p[1:100,:], n_target[1:100],1,2)


# Conclusion 1: realistic RHs are attainable for environmental transmission scenarios
@time βEHμE_mat = get_perc_target(lower_bin,0.05, p, n_target, 11, 15)
@time βEHμH_mat = get_perc_target(lower_bin,0.05, p, n_target, 11, 13)
@time βEHΛH_mat = get_perc_target(lower_bin,0.05, p, n_target, 11, 1)

#plotting!
using ORCA
plotlyjs()
p1 = heatmap(lower_bin, lower_bin, βEHμE_mat, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "μE",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p1.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuEheat.svg")
p2 = heatmap(lower_bin, lower_bin, βEHμH_mat, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "μH",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p2.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuHheat.svg")
p3 = heatmap(lower_bin, lower_bin, βEHΛH_mat, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "ΛH",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p3.o,
        "M:/Project folders/Model env compartment/Plots/bEHLHheat.svg")

#Put into CSVs so I don't have to run again unless I want to
CSV.write("M:/Project folders/Model env compartment/bEHmE.csv", Tables.table(βEHμE_mat))
CSV.write("M:/Project folders/Model env compartment/bEHmH.csv", Tables.table(βEHμH_mat))
CSV.write("M:/Project folders/Model env compartment/bEHLH.csv", Tables.table(βEHΛH_mat))


#Conclusions 2: impact of reducing LA is low for parameter combinations of interest
#I think I should do this for human transmission scenario and then the animal transmission scenario - best and worst chance to have an impact

#human dominated transmission scenario parameters
σ = 1.
MPΛH = 0.1
p_initial[:,1] .= rand(LogNormal(log(MPΛH)+σ, σ), N)
MPβEH = 0.001
σ2 = 2.
p_initial[:,11].= rand(LogNormal(log(MPβEH)+σ2,σ2), N)
MPμE = 0.2
p_initial[:,15] .= rand(LogNormal(log(MPμE)+σ, σ), N)
MPμH = 0.1
p_initial[:,13] .= rand(LogNormal(log(MPμH)+σ, σ), N)

histogram(p_initial[:,[1,11,13,15]],
          xticks = range(0, 1.5; step =0.1),
          xlims = (0.0, 1.5),
          label = ["ΛH" "βEH" "μH" "μE"])

#Rest of the parameters, on tha basis of the
#ΛA and μA
p_initial[:,[2,14]] .= 0.1
#γH, γA, βAA, βAH, βAE, βEA
p_initial[:,[3,4,6,8,9,10]] .= 0.001
#βHH, βHA, βHE
p_initial[:,[5,7,12]] .= 0.2

p = keep_ps(p_initial)

#run with positive lambda value
dat = model_run(p, unboundeds)
#run with 0 lambda value
p[:,2] .= 0
dat_int = model_run(p, unboundeds)

impact = zeros(size(p)[1])
@time for i in 1:size(p)[1]
    if !(dat[i]==0.)
        if dat[i] > dat_int[i]
            impact[i] = (1 - dat_int[i]/dat[i])
        end
    end
end

n_lowimpact = zeros(size(p)[1])
@time for i in 1:size(p)[1]
    if impact[i] < 0.02
        n_lowimpact[i] = 1
    end
end

@time βEHμE_low_impact = get_perc_target(lower_bin,0.05, p, n_lowimpact, 11, 15)
@time βEHμH_low_impact = get_perc_target(lower_bin,0.05, p, n_lowimpact, 11, 13)
@time βEHΛH_low_impact = get_perc_target(lower_bin,0.05, p, n_lowimpact, 11, 1)

p4 = heatmap(lower_bin, lower_bin, βEHμE_low_impact, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "μE",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p4.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuEheat2.svg")
p5 = heatmap(lower_bin, lower_bin, βEHμH_low_impact, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "μH",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p5.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuHheat2.svg")
p6 = heatmap(lower_bin, lower_bin, βEHΛH_low_impact, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "ΛH",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p6.o,
        "M:/Project folders/Model env compartment/Plots/bEHLHheat2.svg")
#################################
#balanced transmission scenario parameters
σ = 1.
MPΛH = 0.1
p_initial[:,1] .= rand(LogNormal(log(MPΛH)+σ, σ), N)
MPβEH = 0.01
σ2 = 2.
p_initial[:,11].= rand(LogNormal(log(MPβEH)+σ2,σ2), N)
MPμE = 0.2
p_initial[:,15] .= rand(LogNormal(log(MPμE)+σ, σ), N)
MPμH = 0.1
p_initial[:,13] .= rand(LogNormal(log(MPμH)+σ, σ), N)

histogram(p_initial[:,[1,11,13,15]],
          xticks = range(0, 1.5; step =0.1),
          xlims = (0.0, 1.5),
          label = ["ΛH" "βEH" "μH" "μE"])

#Rest of the parameters, on tha basis of the
# βEA
p_initial[:,[11]] .= 0.01
#γH, γA, βAH,
p_initial[:,[3,4,7]] .= 0.001
#ΛA, βHH, βAA, βAH, βAE, βHE, μA
p_initial[:,[2,5,6,8,9,12,14]] .= 0.1

p = keep_ps(p_initial)

#run with positive lambda value
dat = model_run(p, unboundeds)
#run with 0 lambda value
p[:,2] .= 0
dat_int = model_run(p, unboundeds)

impact = zeros(size(p)[1])
@time for i in 1:size(p)[1]
    if !(dat[i]==0)
        if dat[i] > dat_int[i]
            impact[i] = (1 - dat_int[i]/dat[i])
        end
    end
end

n_lowimpact = zeros(size(p)[1])
@time for i in 1:size(p)[1]
    if impact[i] < 0.02
        n_lowimpact[i] = 1
    end
end

@time βEHμE_low_impact = get_perc_target(lower_bin,0.05, p, n_lowimpact, 11, 15)
@time βEHμH_low_impact = get_perc_target(lower_bin,0.05, p, n_lowimpact, 11, 13)
@time βEHΛH_low_impact = get_perc_target(lower_bin,0.05, p, n_lowimpact, 11, 1)

p5 = heatmap(lower_bin, lower_bin, βEHμE_low_impact, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "μE",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p4.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuEheat2.svg")
p6 = heatmap(lower_bin, lower_bin, βEHμH_low_impact, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "μH",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p5.o,
        "M:/Project folders/Model env compartment/Plots/bEHmuHheat2.svg")
p7 = heatmap(lower_bin, lower_bin, βEHΛH_low_impact, fillcolor = :fire,
        xlabel = "βEH",
        ylabel = "ΛH",
        xticks = range(0.,stop =1.5, length = 16),
        yticks = range(0.,stop =1.5, length = 16))
PlotlyJS.savefig(p6.o,
        "M:/Project folders/Model env compartment/Plots/bEHLHheat2.svg")


#Final conclusion - increasing βEH in many cases reduces the impacts of increased LH
