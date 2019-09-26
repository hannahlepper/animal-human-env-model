using DifferentialEquations
using Distributions

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
dat = zeros(size(p)[1], 16)
@time for i in 1:size(p)[1]
  dat[i, 1:15] .= p[i,:]
  prob = ODEProblem(unboundeds, u0, tspan, p[i,:])
  sol = solve(prob)
  dat[i,16] = sol(500)[1]
end

dat
