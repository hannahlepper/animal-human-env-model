using DifferentialEquations
using Plots
using GR

#struct so parameters are easier to read/modify
mutable struct ps_struct{T}
  lH    ::T
  lA    ::T
  lE    ::T
  muH   ::T
  muA   ::T
  muE   ::T
  SH    ::T
  SA    ::T
  NH    ::T
  N_A   ::T

  bHH  ::T
  bAA  ::T

  bAH  ::T
  bHA  ::T

  bEH  ::T
  bHE  ::T

  bAE  ::T
  bEA  ::T
end

#fill with values matching B + W 2017 low impact scenario
p = ps_struct{Float16}(
  0.1, 0.1, 0.0,
  0.1, 0.1, 0.0,
  000.0, 000.0,
  1000.0, 1000.0,
  0.1, 0.001,
  0.1, 0.1,
  0.0, 0.0,
  0.0, 0.0)



(1000 - 600) * (0. + (0.1 * 600) / 1000.0 + (0. * 600) / 1000 + 0.1 * 0) - (0.1 * 600)

# differential equations
function an_hum_env(du, u, p, t)
  du[1] =  ((p.NH - u[1]) * (p.lH + (p.bHH * u[1]) / p.NH + (p.bAH * u[2]) / p.NH + p.bEH * u[3])) - p.muH * u[1]
  du[2] =  ((p.N_A - u[2]) * (p.lA + (p.bAA * u[2]) / p.N_A + (p.bHA * u[1]) / p.N_A + p.bEA * u[3])) - p.muA * u[2]
  du[3] =  p.lE + p.SH * p.bHE * u[1] + p.SA * p.bAE * u[2] - p.muE * u[3]
end

# initial values
u0 = [700.0; 700.0; 00.0]

# min and max time steps
tspan = (0.0, 200.0)

#test
an_hum_env_prob = ODEProblem(an_hum_env, u0, tspan, p)
an_hum_env_sol = solve(an_hum_env_prob)

plot(an_hum_env_sol)
an_hum_env_sol(200)

# get Rh*
function run_mod(p)
  prob = ODEProblem(an_hum_env, u0, tspan, p)
  sol = solve(prob)
  sol(200)[1]
end

# empty for RH*
RH = zeros(10000)

# parameters to test
bAH = repeat(collect(0.01:.01:1), outer = 100)
lA = repeat(collect(0.01:.01:1), inner = 100)

for i in 1:10000
  p.bAH = bAH[i]
  p.lA = lA[i]
  RH[i] = run_mod(p)
end

#3d plot
plot(bAH, lA, RH, st = :surface,
       camera = (30,70))

hcat(RH, bAH, lA)
