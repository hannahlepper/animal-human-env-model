# FAST sensitivity analysis
library(fast) #Package gone
library(sensitivity)
library(readr)
library(ggplot2)
library(RColorBrewer)

R_files <- list.files("/Users/hannah/Documents/fast/R")[-"sysdata.rda"]
R_files <- sapply(R_files, function(f) paste(c("/Users/hannah/Documents/fast/R/", f), collapse = ""))
for(f in R_files) source(f)
load("/Users/hannah/Documents/fast/R/sysdata.rda")

#Generate parameter values
p_value_names <- c("LamH", "LamA", "gamH", "gamA",
                   "betHH", "betAA",
                   "betAH", "betEH",
                   "betHA", "betEA",
                   "betHE", "betAE",
                   "muH", "muA", "muE")
para_model = fast_parameters(minimum = c(rep(0, 14), 0.000001),
                             maximum = rep(1, 15),
                             factor = 9,
                             names = p_value_names)
write.csv(para_model, "M:/New projects/Model env compartment/paramodel.csv")

#for wider environmental range
para_model = fast_parameters(minimum = c(rep(0, 14), 0.000001),
                             maximum = c(1, 1, 10, 10, 1, 1, 1, 10, 1, 10, 10, 10, 1, 1, 10),
                             factor = 9,
                             names = p_value_names)
write.csv(para_model, "M:/New projects/Model env compartment/paramodel.csv")

#Equilibrium for each parameter combination in Mathematica
eq_analysis_bounded <- unlist(read_csv("data/outcomeequilibriumBounded.csv", 
                               col_names = FALSE))
eq_analysis_unbounded <- unlist(read_csv("data/outcomeequilibriumUnbounded.csv", 
                               col_names = FALSE))

#sensitivity analysis
sens_bounded <- sensitivity(x = eq_analysis_bounded, 
                    numberf = 15, 
                    make.plot = F, 
                    names = p_value_names)
df.equilibrium_bounded <- data.frame(parameter = factor(p_value_names, ordered = TRUE), 
                             value = sens_bounded)
sens_unbounded <- sensitivity(x = eq_analysis_unbounded, 
                            numberf = 15, 
                            make.plot = F, 
                            names = p_value_names)
df.equilibrium_unbounded <- data.frame(parameter = factor(p_value_names), 
                                     value = sens_unbounded)

p1 <- ggplot(df.equilibrium_bounded, aes(parameter, value))
p2 <- ggplot(df.equilibrium_unbounded, aes(parameter, value))


p2 <- p2 + geom_bar(stat="identity", fill=brewer.pal(3,"Set1")[2]) + 
  scale_x_discrete("Parameter", waiver(),  
                   c(expression(paste(beta[AA])),
                     expression(paste(beta[AE])),
                     expression(paste(beta[AH])),
                     expression(paste(beta[EA])),
                     expression(paste(beta[EH])),
                     expression(paste(beta[HA])),
                     expression(paste(beta[HE])),
                     expression(paste(beta[HH])),
                     expression(paste(gamma[A])),
                     expression(paste(gamma[H])),
                     expression(paste(Lambda[A])),
                     expression(paste(Lambda[H])),
                     expression(paste(mu[A])),
                     expression(paste(mu[E])),
                     expression(paste(mu[H])))) +
  ylab("Partial Variance") + 
    theme_bw() + 
    theme(axis.text=element_text(size = 16, colour = "black"), axis.title=element_text(size=20)) +
  labs(title = "Bounded model FAST")

png("plots/Fig1C.A.png", width = 16, height = 10, units = "cm", res = 300)
p1
dev.off()

png("plots/Fig1C.B.png", width = 16, height = 10, units = "cm", res = 300)
p2
dev.off()
