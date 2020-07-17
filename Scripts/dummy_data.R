#impact data
#rows for randomly drawn parameter set
#columns for transmission scenarios
#elements impact, 0 - 1.0

N <- 10000 #Number of trials in experiments
d_impact <- lapply(1:5, function(ts) rbinom(N, 50, 0.1)/50)


#Parameter data
#List of experiments
#Each experiment a list of transmission scenarios
#Each transmission scenario a df with
# - rows for parameter set
# - columns for the parameter (16 I think)
ts <- c("B", "Bd", "H", "E", "A")
exp1 <- setNames(lapply(1:5, function(ts) {
    data.frame(matrix(rbeta(N*16, 1, 1), ncol = 16))
}), ts)
exp2 <- setNames(lapply(1:5, function(ts) {
    data.frame(matrix(rbeta(N*16, 1, 1), ncol = 16))
}), ts)
Pv <- list(exp1, exp2, exp1, exp2, exp1, exp2)

#bins
lower_bin = seq(0., 1., 0.05)
bin_N = length(lower_bin)