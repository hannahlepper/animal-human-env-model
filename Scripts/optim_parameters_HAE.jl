#Get transmission scenario parameters
using Optim
using DifferentialEquations


function mod_bal(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, bo, β = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + β*RE) - μA*RA
    du[3] = (1 - bo * RE) * (γH*ΛH + γA*ΛA + β*RA + β*RH) - μE*RE
end

function mod_baslowHA(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, βha, β, βE, bo = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + βha*RH + βE*RE) - μA*RA
    du[3] = (1 - bo * RE) * (γH*ΛH + γA*ΛA + β*RA + β*RH) - μE*RE
end

function mod_bashighHA(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, β, βE, bo = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + βE*RE) - μA*RA
    du[3] = (1 - bo * RE) * (γH*ΛH + γA*ΛA + β*RA + β*RH) - μE*RE
end

function mod_hum(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, μH, μA, μE, β, βH, bo = p
    du[1] = (1 - RH) * (ΛH + βH*RH + β*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + βH*RH + β*RE) - μA*RA
    du[3] = (1 - bo * RE) * (γH*ΛH + γA*ΛA + β*RA + βH*RH) - μE*RE
end

function mod_anim(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA,  μH, μA, μE, β, βA, bo = p
    du[1] = (1 - RH) * (ΛH + β*RH + βA*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βA*RA + β*RH + β*RE) - μA*RA
    du[3] = (1 - bo * RE) * (γH*ΛH + γA*ΛA + βA*RA + β*RH) - μE*RE
end

function mod_env(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA,  μH, μA, μE, β, βE, bo = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + βE*RE) - μA*RA
    du[3] = (1 - bo * RE) * (γH*ΛH + γA*ΛA + βE*RA + βE*RH) - μE*RE
end

function dist_targ(model_result) #gets model solution and compares to target - function to be minimised
    abs(0.71 - model_result)
end

function dist_bd(b, bo) # solves balanced model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, bo, b]
    prob = ODEProblem(mod_bal, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_h(b, bo) # solves human model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, b, bo]
    prob = ODEProblem(mod_hum, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_a(b, bo) # solves animal model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, b, bo]
    prob = ODEProblem(mod_anim, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_e(b, bo) # solves env model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, b, bo]
    prob = ODEProblem(mod_env, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_bhighHA(b, bo) #solves baseline 1 model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2,  0.1, b, bo]
    prob = ODEProblem(mod_bashighHA, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end

function dist_blowHA(b, bo) #solves baseline 2 model
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, 1000.)
    p = [0.1, 0.1, 0.001, 0.001, 0.1, 0.1, 0.2, 0.001, 0.1, b, bo]
    prob = ODEProblem(mod_baslowHA, u0, tspan, p)
    sol = solve(prob)
    dist_targ(sol(1000)[1])
end



#optimise for the different transmission scenarios
b_res = map(bo -> optimize(b -> dist_bd(b, bo), 0.0, 1.0).minimizer, 0:1) #0.07432092

blLHA_res = map(bo -> optimize(b -> dist_blowHA(b, bo), 0., 1.).minimizer, 0:1) #0.01521419
blHHA_res = map(bo -> optimize(b -> dist_bhighHA(b, bo), 0., 1.).minimizer, 0:1) #0.003976815

e_res = map(bo -> optimize(b -> dist_e(b, bo), 0., 1.).minimizer, 0:1) #0.1420501
h_res = map(bo -> optimize(b -> dist_h(b, bo), 0., 1.).minimizer, 0:1) #0.2019663
a_res = map(bo -> optimize(b -> dist_a(b, bo), 0., 1.).minimizer, 0:1) #0.2019663
#these agree with plots I have made in Mathematica - happy with the results.
