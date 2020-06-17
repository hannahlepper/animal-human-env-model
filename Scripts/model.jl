using DifferentialEquations

# model function
function unboundeds(du, u, p, t)
    RH, RA, RE = u
    ΛH, ΛA, γH, γA, βHH, βAA, βHA, βAH, βAE, βEA, βEH, βHE, μH, μA, μE = p
    du[1] = (1 - RH) * (ΛH + βHH*RH + βAH*RA + βEH*RE) - μH*RH
    du[2] = (1 - RA) * (ΛA + βAA*RA + βHA*RH + βEA*RE) - μA*RA
    du[3] = γH*ΛH + γA*ΛA + βAE*RA + βHE*RH - μE*RE
end

#plot each transmission scenario
u0 = [0.0; 0.0; 0.0]
tspan = (0.0, 1000.)

#now numerically solve for each parameter set
#function for running model
function model_run(p, mod)
    dat = zeros(size(p)[1],2)
    u0 = [0.0;0.0;0.0]
    tspan = (0., 500.)

    #run 1
    @time for i in 1:size(p)[1]
      prob = ODEProblem(mod, u0, tspan, p[i,:])
      sol = solve(prob)
      dat[i, :] .= Array{Float64,1}(map(n -> sol(n)[1], [400,500]))
    end

    #run 2 - rerun for those possibly not at equilibrium
    not_eqlm = [abs(dat[i,2] - dat[i,1]) > 0.0000001 for i in 1:size(p)[1]]
    rerun = findall(not_eqlm)
    tspan = (0.,10000.)
    @time for i in rerun
        prob = ODEProblem(mod, u0, tspan, p[i,:])
        sol = solve(prob)
        dat[i, :] .= Array{Float64,1}(map(n -> sol(n)[1], [1900,2000]))
    end
    return dat
end


#for testing
p_B =  [0.1, 0.1, 0.001, 0.001, 0.1,        0.1,        0.001,      0.1,        0.1,        0.01,       0.01,       0.1,        0.1, 0.1, 0.2]
@time solve(ODEProblem(unboundeds, u0, tspan, p_B))
