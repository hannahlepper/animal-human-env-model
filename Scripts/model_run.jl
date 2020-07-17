#Get parameters
include("Scripts/parameters.jl")

#Get other processes going
using Distributed
addprocs(11)
nprocs() #check it's working

#Get the model (everyone needs to know)
@time @everywhere include("Scripts/model.jl")

#Generate data
@time dat = pmap(x -> model_run(x, unboundeds), P)
@time dat_orig = pmap(x -> model_run(x, unboundeds), [p_orig, p_orig_int])
#rows for transmission scenarios, columns for experiments

#SAVE/LOAD FILES...
using JLD2
@save "D:/results_workspace.jld2" dat
@save "D:/results_workspace_P.jld2" P
@save "D:/results_workspace_Pv.jld2" Pv
