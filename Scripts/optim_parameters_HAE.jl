#Get transmission scenario parameters
using Optim
using DifferentialEquations


function unboundeds_bal(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, β = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + β*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + β*RA + β*RH - μE*RE
end

function unboundeds_baslowHA(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, βha, β, βE = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + βha*RH + βE*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + β*RA + β*RH - μE*RE
end
function unboundeds_bashighHA(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, β, βE = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + βE*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + β*RA + β*RH - μE*RE
end

function unboundeds_hum(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, β, βH = p
    du[1] = (1 - RH) * (ΛH + βH*RH + β*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + βH*RH + β*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + β*RA + βH*RH - μE*RE
end
function unboundeds_anim(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA,  μH, μA, μE, β, βA = p
    du[1] = (1 - RH) * (ΛH + β*RH + βA*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βA*RA + β*RH + β*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βA*RA + β*RH - μE*RE
end
function unboundeds_env(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA,  μH, μA, μE, β, βE = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + βE*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βE*RA + βE*RH - μE*RE
end

function dist_targ(model_result) #gets model solution and compares to target - function to be minimised
    abs(0.71 - model_result)
end

function dist_bd(b) # solves balanced model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, b]
    prob = ODEProblem(unboundeds_bal, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_h(b) # solves human model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, b]
    prob = ODEProblem(unboundeds_hum, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_a(b) # solves animal model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, b]
    prob = ODEProblem(unboundeds_anim, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_e(b) # solves env model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, b]
    prob = ODEProblem(unboundeds_env, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_bhighHA(b) #solves baseline 1 model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2,  0.1, b]
    prob = ODEProblem(unboundeds_bashighHA, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_blowHA(b) #solves baseline 2 model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, 0.1, b]
    prob = ODEProblem(unboundeds_baslowHA, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end



#optimise for the different transmission scenarios
b_res = optimize(dist_bd, 0., 1.) #0.07432092
blLHA_res = optimize(dist_blowHA, 0., 1.) #0.01521419
blHHA_res = optimize(dist_bhighHA, 0., 1.) #0.003976815
e_res = optimize(dist_e, 0., 1.) #0.1420501
h_res = optimize(dist_h, 0., 1.) #0.2019663
a_res = optimize(dist_a, 0., 1.) #0.2019663
#these agree with plots I have made in Mathematica - happy with the results.
