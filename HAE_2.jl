using DifferentialEquations
using Distributions
using Plots
plotly()
using CSV
using Tables


# model function
function unboundeds(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, βHH, βAA, βHA, βAH, βAE, βEA, βEH, βHE, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + βHH*RH + βAH*RA + βEH*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βAA*RA + βHA*RH + βEA*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βAE*RA + βHE*RH - μE*RE
end

#sample paramter space
N = 1000000
σ = 1.
#parameters of interest are bEH, LH, muE, muH.
p_initial = zeros(N, 15)
#Specify most probably value of LamH, then get lognormal distribution
MPΛH = 0.1
p_initial[:,1] .= rand(LogNormal(log(MPΛH)+σ, σ), N)

MPβEH = 0.14
p_initial[:,11].= rand(LogNormal(log(MPβEH)+σ,σ), N)

MPμE = 0.2
p_initial[:,15] .= rand(LogNormal(log(MPμE)+σ, σ), N)

MPμH = 0.1
p_initial[:,13] .= rand(LogNormal(log(MPμH)+σ, σ), N)

#Rest of the parameters, on tha basis of the
#ΛA and μA
p_initial[:,[2,14]] .= 0.1
#γH, γA, βHH, βAA, βHA, βAH
p_initial[:,[3,4,5,6,7,8]] .= 0.001
#βEA, βAE, βHE
p_initial[:,[9,10,12]] .= 0.14

#get rid of negative numbers, or values over 1? leave this out for now
keep=[]
@time for i in 1:N
    if !(any(p_initial[i,:] .< 0.))
        if !(any(p_initial[i,:] .> 1.5))
            push!(keep, i)
        end
    end
end

p = p_initial[keep, :]

#now numerically solve for each parameter set
dat = zeros(size(p)[1])
@time for i in 1:size(p)[1]
  prob = ODEProblem(unboundeds, u0, tspan, p[i,:])
  sol = solve(prob)
  dat[i] = sol(500)[1]
end

#Find number of runs where 0.65 < RH < 0.75
n_target = zeros(size(p)[1])
for i in 1:size(p)[1]
    if dat[i] < 0.75 && dat[i] < 0.65
        n_target[i] = 1
    else
        n_target[i] = 0
    end
end

#Bins for parameters
maximum(p[:,11])
lower_bin = [0.:0.01:1.;]

#get βEH and μE combinations
βEHμE = zeros(10201,2)
βEHμE[:,1] = repeat(lower_bin; outer=101)
βEHμE[:,2] = repeat(lower_bin; inner=101)

#get indexes of p rows where this is true and calculate %
nd

function get_perc_target(bins, binsize, p, dat, p1, p2)
    dims_out = size(bins)[1]
    dims_in = size(dat)[1]
    perc_target = zeros(dims_out, dims_out)
    for i in 1:dims_out
        for j in 1:dims_out
            sum_success=0
            n_sets=0
            for k in 1:dims_in
                if p[k,p1] > lower_bin[i] && p[k,p1] <  lower_bin[i] + binsize #bEH, starting low
                    if p[k,p2] > lower_bin[j] && p[k,p2] < lower_bin[j] + binsize #muE, startin low
                        n_sets += 1
                        sum_success += n_target[k]
                    end
                end
            end
            perc_target[j,i] = ifelse(sum_success > 0 && n_sets > 0, sum_success/n_sets, 0)
        end
    end
    return perc_target
end
@time βEHμE_mat = get_perc_target(lower_bin,0.01,p, n_target, 11, 15)
@time βEHμH_mat = get_perc_target(lower_bin,0.01,p, n_target, 11, 13)
@time βEHΛH_mat = get_perc_target(lower_bin,0.01,p, n_target, 11, 1)

#plotting!
heatmap(lower_bin, lower_bin, βEHμE_mat)
heatmap(lower_bin, lower_bin, βEHΛH_mat)
heatmap(lower_bin, lower_bin, βEHΛH_mat)

#Put into CSVs so I don't have to run again unless I want to
CSV.write("M:/Project folders/Model env compartment/bEHmE.csv", Tables.table(βEHμE_mat))
CSV.write("M:/Project folders/Model env compartment/bEHmH.csv", Tables.table(βEHμH_mat))
CSV.write("M:/Project folders/Model env compartment/bEHLH.csv", Tables.table(βEHΛH_mat))
