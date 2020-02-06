#Get parameters
include("parameters.jl")

#Get other processes going
using Distributed
addprocs(11)
nprocs() #check it's working

#Get the model (everyone needs to know)
@time @everywhere include("model.jl")

#Generate data
@time dat = pmap(x -> model_run(x, unboundeds), P)
#rows for transmission scenarios, columns for experiments

#SAVE/LOAD FILES...
#using JLD2
#@save "D:/results_workspace.jld2" dat
#@save "D:/results_workspace_P.jld2" P
dat[1,1][:,2]
dat[1,2][:,2]

P[1,1]
P[2,1]
