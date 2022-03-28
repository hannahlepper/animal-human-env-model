
#Get other processes going
using Distributed
addprocs(2)
nprocs() #check it's working

@everywhere (using Pkg; Pkg.activate("."))

#Get the model (everyone needs to know)
using DifferentialEquations
@everywhere using DifferentialEquations
@time @everywhere include("model.jl") 
@time @everywhere include("parameters.jl") 
#sometimes doesn't work - can try closing the workers (rmprocs(2:12)) and reopening them.
#alternatively run on worker 1 first and then the rest

#Generate data
#Unbounded model
@time dat_B = pmap(x -> model_run(get_params(1, x, 2000000, 
        pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:2)

using JLD2
@save "D:/results_workspace_B_1.jld2" dat_B
dat_B = Nothing

@time dat_Bd = pmap(x -> model_run(get_params(2, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:2)
@save "D:/results_workspace_Bd_1.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(3, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:2)
@save "D:/results_workspace_H_1.jld2" dat_H
dat_H = Nothing

@time dat_A = pmap(x -> model_run(get_params(4, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:2)
@save "D:/results_workspace_A_1.jld2" dat_A
dat_A = Nothing

@time dat_E = pmap(x -> model_run(get_params(5, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:2)
@save "D:/results_workspace_E_1.jld2" dat_E
dat_E = Nothing


#Bounded model
@time dat_B = pmap(x -> model_run(get_params(1, x, 2000000, 
        pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:2)
@save "D:/results_workspace_B_2.jld2" dat_B
dat_B = Nothing

@time dat_Bd = pmap(x -> model_run(get_params(2, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:2)
@save "D:/results_workspace_Bd_2.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(3, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:2)
@save "D:/results_workspace_H_2.jld2" dat_H
dat_H = Nothing

@time dat_A = pmap(x -> model_run(get_params(4, x, 2000000, 
    pf2, p_uncertainty2, bEH_experimentsnbounded, LA_experiments), boundeds), 1:2)
@save "D:/results_workspace_A_2.jld2" dat_A
dat_A = Nothing

@time dat_E = pmap(x -> model_run(get_params(5, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:2)
@save "D:/results_workspace_E_2.jld2" dat_E
dat_E = Nothing




@time dat_orig = pmap(x -> model_run(x, unboundeds), [p_orig, p_orig_int])

@load "/mnt/d/results_workspace_P_bo.jld2" P_bo
@time dat_bo = pmap(x -> model_run(x, boundeds), P_bo)
@save "D:/results_workspace_dat_bo.jld2" dat_bo
P_bo = Nothing 


@time dat_orig_bo = pmap(x -> model_run(x, boundeds), [p_orig_bo, p_orig_int_bo])
#rows for transmission scenarios, columns for experiments

#extra data for bounded model results
@time dat_bounded = pmap(x -> model_run(x, unboundeds), [p7])

rmprocs(13:23)
#SAVE/LOAD FILES...
using JLD2
@save "D:/results_workspace.jld2" dat
@save "D:/results_workspace_P.jld2" P
@save "D:/results_workspace_Pv.jld2" Pv
@save "D:/results_workspace_datorig.jld2" dat_orig
