#Get transmission scenario parameters
import("M:/GitHub\\/animal-human-env-model\\/model_parameters_HAE_2.jl")

using Optim


function unboundeds_bal(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, β, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + β*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + β*RA + β*RH - μE*RE
end
function unboundeds_hum(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, β, βH, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + βH*RH + β*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + βH*RH + β*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + β*RA + βH*RH - μE*RE
end
function unboundeds_anim(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, β, βA, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + β*RH + βA*RA + β*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βA*RA + β*RH + β*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βA*RA + β*RH - μE*RE
end
function unboundeds_env(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, β, βE, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + β*RH + β*RA + βE*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + β*RA + β*RH + βE*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βE*RA + βE*RH - μE*RE
end

function solve_mod(b) # selects (based on parameter input) and solves model
    

function diff(b) #gets model solution and compares to target - function to be minimised
