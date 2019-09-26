using DifferentialEquations
using Plots
using PlotlyJS

#struct so parameters are easier to read/modify
mutable struct ps{T}
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

#fill with values matching B + W 2017 low impact scenario
p = ps{Float16}(
  0.1, 0.1,
  0.1, 0.1,
  0.1, 0.1,
  0.1, 0.01,
  0.1, 0.01,
  0.1, 0.1,
  0.1, 0.1, 0.2)

# differential equations
#function an_hum_env(du, u, p, t)
#  du[1] =  ((p.NH - u[1]) * (p.lH + (p.bHH * u[1]) / p.NH + (p.bAH * u[2]) / p.NH + p.bEH * u[3])) - p.muH * u[1]
#  du[2] =  ((p.N_A - u[2]) * (p.lA + (p.bAA * u[2]) / p.N_A + (p.bHA * u[1]) / p.N_A + p.bEA * u[3])) - p.muA * u[2]
#  du[3] =  p.lE + p.SH * p.bHE * u[1] + p.SA * p.bAE * u[2] - p.muE * u[3]
#end

function an_hum_env(du, u, p, t)
  du[1] =  (1 - u[1]) * ((p.lH + (p.bHH * u[1])  + (p.bAH * u[2])  + p.bEH * u[3])) - p.muH * u[1]
  du[2] =  (1 - u[2]) * ((p.lA + (p.bAA * u[2])  + (p.bHA * u[1])  + p.bEA * u[3])) - p.muA * u[2]
  du[3] =  p.gA * p.lA + p.gH * p.lH + p.bHE * u[1] + p.bAE * u[2] - p.muE * u[3]
end

# initial values
u0 = [.000; .000; 00.0]

# min and max time steps
tspan = (0.0, 1000.0)

#test
an_hum_env_prob = ODEProblem(an_hum_env, u0, tspan, p)
an_hum_env_sol = solve(an_hum_env_prob)
plot(an_hum_env_sol)

styles = filter((s->begin
                s in Plots.supported_styles()
            end), [:solid, :dash, :dot, :dashdot, :dashdotdot])
styles = reshape(styles, 1, length(styles))

p.bHA = 0.1
an_hum_env_prob = ODEProblem(an_hum_env, u0, tspan, p)
an_hum_env_sol = solve(an_hum_env_prob)
p1 = Plots.plot(an_hum_env_sol[[1,2],:]', line = (2, styles),
           label = [:RH, :RA],
           title = "No environment, low impact",
           yaxis = ("", (0, 1.0), (0:0.1:1.0)))

p.bHA = 0.001
an_hum_env_sol = solve(an_hum_env_prob)
p2 = Plots.plot(an_hum_env_sol[[1,2],:]', line = (2, styles),
           label = [:RH, :RA],
           title = "No environment, high impact",
           yaxis = ("", (0, 1.0), (0:0.1:1.0)))


p.bAE = 0.1
p.bHE = 0.1
p.bEH = 0.01
p.bEA = 0.01
p.muE = 0.1
p.lE = 0.1
p.bHA = 0.001
an_hum_env_prob = ODEProblem(an_hum_env, u0, tspan, p)
an_hum_env_sol = solve(an_hum_env_prob)
p3 = Plots.plot(an_hum_env_sol[[1,2],:]', line = (2, styles),
           label = [:RH, :RA],
           title = "Less environment, high impact",
           yaxis = ("", (0, 1.0), (0:0.1:1.0)))

p.bHA = 0.1
an_hum_env_sol = solve(an_hum_env_prob)
p4 = Plots.plot(an_hum_env_sol[[1,2],:]', line = (2, styles),
           label = [:RH, :RA],
           title = "Less environment, low impact",
           yaxis = ("", (0, 1.0), (0:0.1:1.0)))

p1
p2
p3
p4
p5
p6

p.bAE = 0.1
p.bHE = 0.1
p.bEH = 0.1
p.bEA = 0.1
p.muE = 0.1
p.lE = 0.1

p.bHA = 0.1
an_hum_env_sol = solve(an_hum_env_prob)
p5 = Plots.plot(an_hum_env_sol[[1,2],:]', line = (2, styles),
           label = [:RH, :RA],
           title = "High environment, low impact",
           yaxis = ("", (0, 1.0), (0:0.1:1.0)))

p.bHA = 0.001
an_hum_env_sol = solve(an_hum_env_prob)
p6 = Plots.plot(an_hum_env_sol[[1,2],:]', line = (2, styles),
          label = [:RH, :RA],
          title = "High environment, high impact",
          yaxis = ("", (0, 1.0), (0:0.1:1.0)))

# get Rh*
function run_mod(p)
  prob = ODEProblem(an_hum_env, u0, tspan, p)
  sol = solve(prob, alg_hints = [:stiff])
  sol(200)[1]
end

# function to get the simulated variables
function sens_run()
  RH = zeros(121)
  bAH = repeat(collect(0.0:.1:1), inner = 11)
  lA = repeat(collect(0.0:.1:1), outer = 11)
  for i in 1:121
    p.bAH = bAH[i]
    p.lA = lA[i]
    RH[i] = run_mod(p)
  end
  return RH
end

# function to get 3d plot

function plotlyplot3d(RH, title)
  RH_vec = Vector([RH[i:i+10] for i in collect(1:11:121)])
  trace_RH = PlotlyJS.surface(z=RH_vec,
                              x = collect(0.0:.1:1),
                              y = collect(1.0:-.1:0.0),
                              colorscale="Viridis",
                              autoscale = false,
                              )
  layout_RH = PlotlyJS.Layout(;title = title)
  PlotlyJS.plot(trace_RH, layout_RH)
end

# Scenario 1

p.bAE = 0.0
p.bHE = 0.0
p.bEH = 0.0
p.bEA = 0.0
p.muE = 0.0
p.lE = 0.0

p.bHA = 0.1

RH = sens_run()
p1sc1 = plotlyplot3d(RH, "no environment, low impact")

# Scenario 2

p.bHA = 0.001
RH = sens_run()
p2sc2 = plotlyplot3d(RH, "no environment, high impact")

[p1sc1, p2sc2,
 p3sc1, p4sc2,
 p5sc1, p6sc2]
# Scenario 1; some environment

p.bHA = 0.1
p.bAE = 0.1
p.bHE = 0.1
p.bEA = 0.01
p.bEH = 0.01
p.muE = 0.5
p.lE = 0.01

RH = sens_run()
p3sc1 = plotlyplot3d(RH, "some environment, low impact")

# Scenario 2; some environment

p.bHA = 0.01

RH = sens_run()
p4sc2 = plotlyplot3d(RH, "some environment, high impact")


# Scenario 1; high environment

p.bHA = 0.1
p.bAE = 0.1
p.bHE = 0.1
p.bEA = 0.1
p.bEH = 0.1
p.muE = 0.5
p.lE = 0.01

RH = sens_run()
p5sc1 = plotlyplot3d(RH, "high environment, low impact")


# Scenario 2; high environment

p.bHA = 0.01

RH = sens_run()
p6sc2 = plotlyplot3d(RH, "high environment, high impact")


# FAST sensitivity analysis
using CSV
using DataFrames

ps(xs) = ps{Float64}(xs[1],xs[2],xs[3],xs[4],xs[5],xs[6],xs[7],xs[8],xs[9],xs[10],xs[11],xs[12],xs[13],xs[14],xs[15])

# initial values
u0 = [.000; .000; 00.0]
# min and max time steps
tspan = (0.0, 1000.0)

function run_mod(p)
  prob = ODEProblem(an_hum_env, u0, tspan, p)
  sol = solve(prob, alg_hints = [:stiff])
  sol(1000)[1]
end


par_df = CSV.read("/Users/hannah/Dropbox/Academic/Edinburgh/paramodel.csv";
         header = 1,
         delim = ",",
         types = fill(Float64, 16))
par_mat = Matrix(par_df)
RH = [run_mod(ps(par_mat[xs,2:16])) for xs in 1:9171]
CSV.write("/Users/hannah/Dropbox/Academic/Edinburgh/equilibrium_paramodel.csv", DataFrame(RH = RH),
          delim = ",", writeheader = false)
