
#Get other processes going
using Distributed
addprocs(9)
nprocs() #check it's working

@everywhere (using Pkg; Pkg.activate("."))

#Get the model
using DifferentialEquations
using JLD2
@everywhere using DifferentialEquations
@time @everywhere include("model.jl") 
@time @everywhere include("parameters.jl") 
#sometimes doesn't work - can try closing the workers (rmprocs(2:12)) and reopening them.
#alternatively run on worker 1 first and then the rest

#Generate data
#Unbounded model
N_sets = 1000000
@time dat_Bd = pmap(x -> model_run(get_params(1, x, N_sets, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_revised_Bd_1.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(2, x, N_sets, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_revised_H_1.jld2" dat_H
dat_H = Nothing

@time dat_E = pmap(x -> model_run(get_params(3, x, N_sets, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_revised_E_1.jld2" dat_E
dat_E = Nothing

@time dat_A = pmap(x -> model_run(get_params(4, x, N_sets, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_revised_A_1.jld2" dat_A
dat_A = Nothing


#Bounded model

@time dat_Bd = pmap(x -> model_run(get_params(1, x, N_sets, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_revised_Bd_2.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(2, x, N_sets, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_revised_H_2.jld2" dat_H
dat_H = Nothing

@time dat_E = pmap(x -> model_run(get_params(3, x, N_sets, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_revised_E_2.jld2" dat_E
dat_E = Nothing

@time dat_A = pmap(x -> model_run(get_params(4, x, N_sets, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_revised_A_2.jld2" dat_A
dat_A = Nothing

#Original model

@time dat_Bd = pmap(x -> model_run(get_params(1, x, N_sets, 
    pf3, p_uncertainty3, bEH_experiments_orig, LA_experiments), orig), 1:13)
@save "/mnt/d/results_workspace_revised_Bd_3.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(2, x, N_sets, 
    pf3, p_uncertainty3, bEH_experiments_orig, LA_experiments), orig), 1:13)
@save "/mnt/d/results_workspace_revised_H_3.jld2" dat_H
dat_H = Nothing

@time dat_A = pmap(x -> model_run(get_params(4, x, N_sets, 
    pf3, p_uncertainty3, bEH_experiments_orig, LA_experiments), orig), 1:13)
@save "/mnt/d/results_workspace_revised_A_3.jld2" dat_A
dat_A = Nothing
