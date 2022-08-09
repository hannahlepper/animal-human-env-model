#Get transmission scenario parameters
using Optim
using DifferentialEquations

function dist_targ(model_targets, model_results) #gets model solution and compares to target - function to be minimised
    #sum(abs2, model_results - model_targets)
    (model_results - model_targets)^2
end

function mod(du, u, p, t)
    RH, RA, RE = u
    γH, γA, μH, μA, μE, ΛH, ΛA, βHH, βAA, βAH, βHA, βEH, βEA, βHE, βAE, bound, orig = p
    du[1] = (1 - RH) * (ΛH + βHH*RH + βAH*RA + orig*βEH*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βAA*RA + βHA*RH + orig*βEA*RE) - μA*RA
    du[3] = orig * ((1 - bound * RE) * (γH*ΛH + γA*ΛA + βAE*RA + βHE*RH) - μE*RE)
end

function dist(x, targ, mod, bound, orig, TS, return_sol)
    
    γH, γA, μH, μA, μE, ΛH, ΛA = [0.001, 0.001, 0.118, 0.118*2, 0.29, 0.1, 0.1] #specify constant parameters

    #Specify parameters for different transmission scenarios
    if TS == :balance
        βHH = βAA = βAH = βHA = βEH = βEA = βHE = βAE = x  #All parameters to be the minimizer
    elseif TS == :human
        βHH = βHA = βHE = x #Some parameters to be the minimizer
        βAA = βAH = βEH = βEA = βAE = 0.001 #Others to be set low
    elseif TS == :animal
        βAA = βAH = βAE = x 
        βHH = βHA = βEH = βHE = βEA = 0.001
    elseif TS == :environment
        βEH = βEA = βHE = βAE = x  
        βHH = βAA = βAH = βHA = 0.001
    end

    p = [γH, γA, μH, μA, μE, ΛH, ΛA, βHH, βAA, βAH, βHA, βEH, βEA, βHE, βAE, bound, orig]
    u0 = [0.0; 0.0; 0.0] #initial conditions
    tspan = (0.0, 1000.) #time span
    prob = ODEProblem(mod, u0, tspan, p) #define model as correct type
    sol = solve(prob) #obtain model solution
    if return_sol
        sol = solve(prob) #return solution (for checking)
    else
        dist_targ(targ, sol(1000)[1]) #or return distance between target and solution (for fitting)
    end
end


#optimise for the different transmission scenarios
b_bounded = map(ts -> optimize(b -> dist(b, 0.541, mod, 1, 1, ts, false), 0.0, 1.0).minimizer, [:balance, :human, :animal, :environment]) 
b_unbounded = map(ts -> optimize(b -> dist(b, 0.541, mod, 0, 1, ts, false), 0.0, 1.0).minimizer, [:balance, :human, :animal, :environment]) 
b_orig = map(ts -> optimize(b -> dist(b, 0.541, mod, 0, 0, ts, false), 0.0, 1.0).minimizer, [:balance, :human, :animal, :environment]) 

#             balance    human   animal   environment
#bounded       0.0188   0.0315   0.0510        0.0798
#unbounded     0.0188   0.0315   0.0510        0.0741
#original      0.0200   0.0316   0.0512             -

b_bounded_bal, b_bounded_hum, b_bounded_anim, b_bounded_env = b_bounded
b_unbounded_bal, b_unbounded_hum, b_unbounded_anim, b_unbounded_env = b_unbounded
b_orig_bal, b_orig_hum, b_orig_anim, b_orig_env = b_orig

transmission_scenario = [:balance, :human, :animal, :environment]
map(ts -> dist(b_bounded[ts], 0.541, mod, 1, 1, transmission_scenario[ts], true)(1000), 1:4)
map(ts -> dist(b_unbounded[ts], 0.541, mod, 0, 1, transmission_scenario[ts], true)(1000), 1:4)
map(ts -> dist(b_orig[ts], 0.541, mod, 0, 0, transmission_scenario[ts], true)(1000), 1:4)