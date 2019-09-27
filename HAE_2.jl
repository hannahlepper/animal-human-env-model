using DifferentialEquations
using Distributions
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

#dummy run to show functionality
p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
u0 = [0.0;0.0;0.0]
tspan=(0.0,500.0)
prob = ODEProblem(unboundeds,u0,tspan,p)
@time sol = solve(prob)
sol(500)

#sample paramter space
N = 1000
σ = 1
#parameters of interest are bEH, LH, muE, muH.
p_initial = zeros(1000, 15)
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
@time for i in 1:1000
    if !(any(p_initial[i,:] .< 0.))
        #if !(any(p_initial[i,:] .> 1.))
            push!(keep, i)
        #end
    end
end

p = p_initial[keep, :]

#now solve for each parameter set
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
maximum(p[:,11]) #too big - not interesting this high up in the range
lower_bin = [0.:0.2:1.8;]

#get βEH and μE combinations
βEHμE = zeros(100,2)
βEHμE[:,1] = repeat(lower_bin; outer=10)
βEHμE[:,2] = repeat(lower_bin; inner=10)

#get indexes of p rows where this is true and calculate %
perc_target = []
for i in 1:100
    #go through the parameter combinations first
    indexes_i = []
    for j in 1:size(p)[1]
        #go through each parameter set
        if p[j,11] > βEHμE[i,1] && p[j,11] <  βEHμE[i,1] + 0.2
            if p[j,15] > βEHμE[i,2] && p[j,15] < βEHμE[i,2] + 0.2
                push!(indexes_i, j)
            end
        end
    end
    if any(n_target[indexes_i] .> 0)
        sum_success = sum(n_target[indexes_i])
        perc = sum_success/size(indexes_i)[1]
        push!(perc_target, perc)
    else
        push!(perc_target, 0)
    end
end

db = Tables.table(hcat(βEHμE, perc_target))

CSV.write("betEHmuEresults.csv",db)
