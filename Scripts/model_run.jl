#Get parameters
include("Scripts/parameters.jl")

#Get other processes going
using Distributed
addprocs(11)
nprocs() #check it's working

#Get the model (everyone needs to know)
@time @everywhere include("Scripts/model.jl") 
#sometimes doesn't work - can try closing the workers (rmprocs(2:12)) and reopening them.
#alternatively run on worker 1 first and then the rest

#Generate data
@time dat = pmap(x -> model_run(x, unboundeds), P)
@time dat_orig = pmap(x -> model_run(x, unboundeds), [p_orig, p_orig_int])
#rows for transmission scenarios, columns for experiments

rmprocs(13:23)
#SAVE/LOAD FILES...
using JLD2
@save "D:/results_workspace.jld2" dat
@save "D:/results_workspace_P.jld2" P
@save "D:/results_workspace_Pv.jld2" Pv
@save "D:/results_workspace_datorig.jld2" dat_orig
