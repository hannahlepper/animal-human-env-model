using DifferentialEquations
using Plots

#struct so parameters are easier to read/modify
mutable struct ps_struct{T}
 lH    ::T
 lA    ::T

 gH    ::T
 gA    ::T

 bHH  ::T
 bAA  ::T

 bAH  ::T
 bEH  ::T

 bHA  ::T
 bEA  ::T

 bHE  ::T
 bAE  ::T

 muH   ::T
 muA   ::T
 muE   ::T
end

#fill with base values
p = ps_strcut{Float16}(
 0.1, 0.1,
 0.1, 0.1,
 0.1, 0.1,
 0.1, 0.01,
 0.1, 0.01,
 0.1, 0.1,
 0.1, 0.1, 0.2)

# differential equations
function an_hum_env(du, u, p, t)
 du[1] =  (1 - u[1]) * (p.lH + p.bHH * u[1] + p.bAH * u[2] + p.bEH * u[3]) - p.muH * u[1]
 du[2] =  (1 - u[2]) * (p.lA + p.bAA * u[2] + p.bHA * u[1] + p.bEA * u[3]) - p.muA * u[2]
 du[3] =  p.gA * p.lA + p.gH * p.lH + p.bHE * u[1] + p.bAE * u[2] - p.muE * u[3]
end

function an_hum_env_b(du, u, p, t)
  du[1] =  ((1 - u[1]) * (p.lH + (p.bHH * u[1]) + (p.bAH * u[2]) + p.bEH * u[3])) - p.muH * u[1]
  du[2] =  ((1 - u[2]) * (p.lA + (p.bAA * u[2]) + (p.bHA * u[1]) + p.bEA * u[3])) - p.muA * u[2]
  du[3] =  ((1 - u[3]) * (p.gA * p.lA + p.gH * p.lH + p.bHE * u[1] + p.bAE * u[2])) - p.muE * u[3]
end

# initial values
u0 = [.0; .0; .0]

# min and max time steps
tspan = (0.0, 200.0)

#test
an_hum_env_prob = ODEProblem(an_hum_env, u0, tspan, p)
an_hum_env_sol = solve(an_hum_env_prob)
plot(an_hum_env_sol)

an_hum_env_b_prob = ODEProblem(an_hum_env_b, u0, tspan, p)
an_hum_env_b_sol = solve(an_hum_env_b_prob)
plot(an_hum_env_b_sol)

# get Rh*
function run_mod(p)
  prob = ODEProblem(an_hum_env, u0, tspan, p)
  sol = solve(prob, alg_hints = [:stiff])
  sol(200)[1]
end

# empty for RH*
RH = zeros(121)

# parameters to test
bAH = repeat(collect(0.0:.1:1), outer = 11)
lA = repeat(collect(0.0:.1:1), inner = 11)

for i in 1:121
  p.bAH = bAH[i]
  p.lA = lA[i]
  RH[i] = run_mod(p)
end

#3d plot

RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
trace_RH = surface(z=RH_vec, colorscale="Viridis")
plot(trace_RH)

# Scenario 2

p.bHA = 0.001

for i in 1:121
  p.bAH = bAH[i]
  p.lA = lA[i]
  RH[i] = run_mod(p)
end

RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
trace_RH = surface(z=RH_vec, colorscale="Viridis")
plot(trace_RH)

# Scenario 1; some environment

p.bHA = 0.1
p.bAE = 0.001
p.bHE = 0.001
p.bEA = 0.001
p.bEH = 0.001
p.muE = 0.5
p.lE = 0.001

for i in 1:121
  p.bAH = bAH[i]
  p.lA = lA[i]
  RH[i] = run_mod(p)
end

RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
trace_RH = surface(z=RH_vec, colorscale="Viridis")
plot(trace_RH)

# Scenario 2; some environment

p.bHA = 0.01
p.bAE = 0.001
p.bHE = 0.001
p.bEA = 0.001
p.bEH = 0.001
p.muE = 0.5
p.lE = 0.001

for i in 1:121
  p.bAH = bAH[i]
  p.lA = lA[i]
  RH[i] = run_mod(p)
end

RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
trace_RH = surface(z=RH_vec, colorscale="Viridis")
plot(trace_RH)

# Scenario 1; high environment

p.bHA = 0.1
p.bAE = 0.1
p.bHE = 0.1
p.bEA = 0.1
p.bEH = 0.1
p.muE = 0.5
p.lE = 0.001

for i in 1:121
 p.bAH = bAH[i]
 p.lA = lA[i]
 RH[i] = run_mod(p)
end

RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
trace_RH = surface(z=RH_vec, colorscale="Viridis")
plot(trace_RH)

# Scenario 2; high environment

p.bHA = 0.01
p.bAE = 0.1
p.bHE = 0.1
p.bEA = 0.1
p.bEH = 0.1
p.muE = 0.5
p.lE = 0.001

for i in 1:121
 p.bAH = bAH[i]
 p.lA = lA[i]
 RH[i] = run_mod(p)
end

RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
trace_RH = surface(z=RH_vec, colorscale="Viridis")
plot(trace_RH)
