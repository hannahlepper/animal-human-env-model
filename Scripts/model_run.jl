
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
@time dat_B = pmap(x -> model_run(get_params(1, x, 2000000, 
        pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_B_1.jld2" dat_B
dat_B = Nothing

@time dat_Bd = pmap(x -> model_run(get_params(2, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_Bd_1.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(3, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_H_1.jld2" dat_H
dat_H = Nothing

@time dat_E = pmap(x -> model_run(get_params(4, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_A_1.jld2" dat_E
dat_E = Nothing

@time dat_A = pmap(x -> model_run(get_params(5, x, 2000000, 
    pf, p_uncertainty, bEH_experiments_unbounded, LA_experiments), unboundeds), 1:13)
@save "/mnt/d/results_workspace_E_1.jld2" dat_A
dat_A = Nothing


#Bounded model
@time dat_B = pmap(x -> model_run(get_params(1, x, 2000000, 
        pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_B_2.jld2" dat_B
dat_B = Nothing

@time dat_Bd = pmap(x -> model_run(get_params(2, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_Bd_2.jld2" dat_Bd
dat_Bd = Nothing

@time dat_H = pmap(x -> model_run(get_params(3, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_H_2.jld2" dat_H
dat_H = Nothing

@time dat_E = pmap(x -> model_run(get_params(4, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_A_2.jld2" dat_E
dat_A = Nothing

@time dat_A = pmap(x -> model_run(get_params(5, x, 2000000, 
    pf2, p_uncertainty2, bEH_experiments_bounded, LA_experiments), boundeds), 1:13)
@save "/mnt/d/results_workspace_E_2.jld2" dat_A
dat_A = Nothing



#Data for original models
#@load "/mnt/d/results_workspace_p_orig.jld2" p_orig
#@load "/mnt/d/results_workspace_p_orig_bo.jld2" p_orig_bo
#@load "/mnt/d/results_workspace_p_orig_int.jld2" p_orig_int
#@load "/mnt/d/results_workspace_p_orig_int_bo.jld2" p_orig_int_bo

@time dat_orig_unbounded = pmap(x -> model_run(x, unboundeds), [p_orig, p_orig_int])
@time dat_orig_bounded = pmap(x -> model_run(x, boundeds), [p_orig_bo, p_orig_int_bo])

@save "/mnt/d/results_workspace_datorig.jld2" dat_orig_unbounded
@save "/mnt/d/results_workspace_datorig_bo.jld2" dat_orig_bounded
