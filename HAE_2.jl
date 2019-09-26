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

#sample paramter space - dummy
p_initial = zeros(1000, 15)
for i in 1:15
    p_initial[:,i] = rand(LogNormal(log(0.1), 0.05), 1000)
end

#get rid of negative numbers
keep=[]
for i in 1:1000
    if !(any(p_initial[i,:] .< 0.))
        if !(any(p_initial[i,:] .> 1.))
            push!(keep, i)
        end
    end
end

p = p_initial[keep, :]

#now solve for each parameter set
dat = zeros(size(p)[1], 16)
for i in 1:size(p)[1]
  dat[i, 1:15] .= p[i,:]
  prob = ODEProblem(unboundeds, u0, tspan, p[i,:])
  sol = solve(prob)
  dat[i,16] = sol(500)[1]
end

dat
