
#Get other processes going
using Distributed
addprocs(9)
nprocs() #check it's working

@everywhere (using Pkg; Pkg.activate("."))

#Get the model (everyone needs to know)
using DifferentialEquations
@everywhere using DifferentialEquations
@time @everywhere include("model.jl") 
#sometimes doesn't work - can try closing the workers (rmprocs(2:12)) and reopening them.
#alternatively run on worker 1 first and then the rest

#Generate data
using JLD2

@load "/mnt/d/results_workspace_P.jld2" P
@time dat = pmap(x -> model_run(x, unboundeds), P)
@save "D:/results_workspace.jld2" dat
p = Nothing

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
