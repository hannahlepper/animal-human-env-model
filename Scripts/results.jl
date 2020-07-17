# using JLD2
#
# @load "D:/results_workspace.jld2" dat
# @load "D:/results_workspace_P.jld2" P


#GET PRESENCE/ABSENCE OF TARGETS REACHED

#Get measures of impact
function get_impact(dat_pre_int, dat_post_int)
    if !(dat_pre_int==0.) #can't have a 0 in denominator
        1 - dat_post_int/dat_pre_int #will get %decrease and %increase here
    else
        1 - dat_post_int/(dat_pre_int + 0.0000001) #will get %decrease and %increase here
    end
end

function map_impacts(dat_pre_col, dat_post_col)
    map(x -> get_impact(dat_pre_col[x], dat_post_col[x]), 1:length(dat_pre_col))
end

#1.fixed LA intervention, varying bEH
impact_1 = [map_impacts(dat[i,1][:,2], dat[i,2][:,2]) for i in 1:5]

#2. varyng LA intervention, varying bEH
impact_2 = [map_impacts(dat[i,3][:,2], dat[i,2][:,2]) for i in 1:5]

#3. varying bEH intervention, fixed LA
impact_3 = [map_impacts(dat[i,1][:,2], dat[i,5][:,2]) for i in 1:5]

#4. varying LA intervention, fixed bEH
impact_4 = [map_impacts(dat[i,4][:,2], dat[i,6][:,2]) for i in 1:5]

#5. varying LA intervention, varying bEH, low bHA
impact_5 = [map_impacts(dat[i,9][:,2], dat[i, 8][:,2]) for i in 1:5]

#6. fixed bEH intervention
impact_6 = [map_impacts(dat[i,7][:,2], dat[i,5][:,2]) for i in 1:5]

#7. fixed LA intervention
impact_7 = [map_impacts(dat[i,7][:,2], dat[i,6][:,2]) for i in 1:5]

#8. both bEH and LA reduced to 0
impact_8 = [map_impacts(dat[i,7][:,2], dat[i, 10][:,2]) for i in 1:5]

#9. impact of original model - high impact scenario
impact_9 = map_impacts(dat_orig[1][:,2], dat_orig[2][:,2])
mean(impact_9)

#Did simulations reach target RH of 0.65 - 0.75?
n_target = hcat([Int.(0.65 .< dat[x,1][:,1] .< 0.75) for x in 1:5])

#Did simulations have a low impact of less than 2%?
n_low_impact = hcat([Int.(impact_1[x] .< 0.02) for x in 1:5])


#Using RCall to get % reaching target
using RCall
R"source('RGetPercAndPlot.R')"

@rput Pv

#Fig 1 B
#Bar plot of RH, RA, and RE values - deterministic results
R = pmap(x -> solve(ODEProblem(unboundeds, u0, tspan, p[x]))(200), 1:5)
@rput R
R"source('Scripts/Fig 1/RHforeachTSplot.R')"

#Fig 2 A and B 
#Compare bEH and LA interventions
@rput Pv impact_3 impact_4
R"source('Scripts/Fig 2/fig2.R')"

#Fig 3 A, B and apendix
#Heatmaps of impact as LA and bEH changes, for two levels of bHA
#Bins for parameters
lower_bin = [0.:0.05:1;]
bin_N = size(lower_bin)[1]
@rput impact_2 impact_5
R"source('Scripts/Fig 3/Fig3A.R')"







# Conclusion 1: realistic RHs are attainable for environmental transmission scenarios
R"""
mean_impact = lapply(c(15,13,1), function(p) {
        lapply(1:5, function(x) {
            get_perc_target(lower_bin, Pv[[1]][[x]][, c(11,p)], n_target[[x]])
        })
})
names(mean_impact) = c("muE", "muH", "LH")

#plotting!
bEH_mean_impact_plots = lapply(c("muE", "muH", "LH"), function(p) {
    lapply(1:5, function(x) {
        plot_heatmap(mean_impact[[p]][[x]], c('bEH', p),
                     paste(c("% Target Achieved, ts = ", c("B", "Bd", "H", "E", "A")[x]), collapse = ""),
        limits = c(0, 0.75))
        })
    })
"""
R"""
#svg('M:/Project folders/Model env compartment/Plots/ptaplot.svg', height=20, width = 20)
wrap_plots(unlist(bEH_mean_impact_plots, recursive = FALSE), ncol = 5)
#dev.off()
"""

#Conclusions 2: impact of reducing LA is low for parameter combinations of interest
R"""
low_impact_prop = lapply(c(15,13,1), function(p) {
        lapply(1:5, function(x) {
            get_perc_target(lower_bin, Pv[[1]][[x]][, c(11,p)], n_low_impact[[x]])
        })
})
names(low_impact_prop) = c("muE", "muH", "LH")

#plotting!
low_impact_prop_plots = lapply(c("muE", "muH", "LH"), function(p) {
    lapply(1:5, function(x) {
        plot_heatmap(low_impact_prop[[p]][[x]], c('bEH', p),
                     paste(c("% low (<2%) impact, ts = ", c("B", "Bd", "H", "E", "A")[x]), collapse = ""),
        limits = c(0, 1))
        })
    })
"""
R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/liplot.svg', height=20, width = 20)
low_impact_plots = wrap_plots(unlist(low_impact_prop_plots, recursive = FALSE), ncol = 5)
#dev.off()
"""


#Conclusion 3: increasing βEH in many cases reduces the impacts of LA interventions
R"""
mean_impact <- lapply(1:5, function(ts) {
    get_mean_var(lower_bin, Pv[[3]][[ts]][,c(11,2)], impact_2[ts])
})
names(mean_impact) <- c("B", "Bd", "H", "E", "A")
mean_impact_plot_list = lapply(1:5, function(ts) {
    plot_heatmap(mean_impact[[ts]], c('bEH', 'LA'),
    paste(c("Mean impact, ts = ", c("B", "Bd", "H", "E", "A")[ts]), collapse = ""),
    limits = c(0., 0.126))
})
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miLAplot.svg', height=20, width = 5)
miLAplot = wrap_plots(mean_impact_plot_list)
#dev.off()
"""
R"""dev.off()"""








#What is the effect of lower bHA? Do we replicate the results of the original study?
R"""
mean_impact <- lapply(1:5, function(ts) {
    get_mean_var(lower_bin, Pv[[9]][[ts]][,c(11,2)], impact_5[ts])
})
names(mean_impact) <- c("B", "Bd", "H", "E", "A")

mean_impact_plot_list = lapply(1:5, function(ts) {
    plot_heatmap(mean_impact[[ts]], c('bEH', 'LA'),
    paste(c("Low bHA, mean impact, ts = ", c("B", "Bd", "H", "E", "A")[ts]), collapse = ""),
    limits = c(0., 0.126))
})
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miLAplot.svg', height=20, width = 5)
miLAplot = wrap_plots(mean_impact_plot_list)
#dev.off()
"""
R"""dev.off()"""

#Compare bEH and LA interventions
R"""
dfs <- adply(1:5, 1, function(ts) {
    data.frame(param_val = c(Pv[[1]][[ts]][1:10000,11], Pv[[4]][[ts]][1:10000,2]),
               param = rep(c("bEH", "LA"), each = 10000),
               impact = c(impact_3[[ts]][1:10000], impact_4[[ts]][1:10000]))
})
head(dfs)
dfs = mutate(dfs, ts = rep(c("B", "Bd", "H", "E", "A"), each = 20000)) %>%
      mutate(., pv_cut = cut(dfs$param_val, c(0, 0.095, 0.105, 0.495, 0.505, 1.0)))
df_low <- dfs[dfs$pv_cut == "(0.095,0.105]",]
df_high <- dfs[dfs$pv_cut == "(0.495,0.505]",]

wrap_plots(ggplot(df_low, aes(ts, impact, col = param)) + geom_boxplot(),
ggplot(df_high, aes(ts, impact, col = param)) + geom_boxplot())
"""

#Get together RH measures for these indices
R"""
RH <- data.frame(ts =  rep(c("B", "Bd", "H", "E", "A"), each = n*4),
                 bEH = rep(c(0.14, 0.01, 0.07432092, 0.0001, 0.0001), each = 4),
                 RH = c(dat_E_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_B_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_Bd_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_A_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_H_7[sample(1:190000,10000, replace = FALSE),2]))

svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/TransmissionScenarioRHVar.svg', height=7, width = 7)
ggplot(RH, aes(y=ts, x=RH, fill = factor(bEH))) +
    geom_density_ridges(scale=1) +
    #geom_violin(outlier.shape = NA, kernel = "rectangular") +
    #geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_vline(xintercept = 0.71) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour= "black"),
          panel.grid = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          strip.background = element_rect(fill = NA, colour = NA),
          axis.text = element_text(size = 15))
dev.off()
"""

R"dev.off()"


#compare impacts of 0.1 -> 0.0over bEH variation across TS
R"""
dfs <- adply(1:5, 1, function(ts) {
    data.frame(bEH_val = c(Pv[[1]][[ts]][,11]),
               impact = c(impact_1[[ts]]))
})
str(dfs)
dfs_orig <- data.frame(X1 = rep(1, length(impact_1[[1]])),
                        bEH_val = rep(0, length(impact_1[[1]])),
                       impact = impact_9)
dfs = rbind(dfs, dfs_orig)
dfs = mutate(dfs, ts = rep(c('B', 'Bd', 'H', 'E', 'A', 'orig'), each = length(impact_1[[1]]))) %>%
      mutate(., pv_cut = cut(dfs$bEH_val, seq(0., 1., 0.01)))
dfs_summary <- ddply(dfs, .(ts, pv_cut), summarise, mean_impact = mean(impact), sd_impact = prop_sd(mean(impact), length(impact)), bEH = mean(bEH_val))

ggplot(dfs_summary, aes(bEH, mean_impact, col = ts)) + geom_point(position = 'jitter') +
    #scale_y_continuous(limits = c(0, 0.3)) +
    geom_errorbar(aes(ymin = mean_impact - sd_impact, ymax = mean_impact + sd_impact)) +
    geom_line()
"""
#Compare only bEH or LA interventions and bEH + LA intervention
R"""
#need impacts 6, 7 and 8 - balanced only
mean_bEH_only = sapply(1:5, function(ts) mean(impact_6[[ts]]))
mean_LA_only = sapply(1:5, function(ts) mean(impact_7[[ts]]))
mean_both = sapply(1:5, function(ts) mean(impact_8[[ts]]))
sd_bEH_only = sapply(1:5, function(ts) sd(impact_6[[ts]]))
sd_LA_only = sapply(1:5, function(ts) sd(impact_7[[ts]]))
sd_both = sapply(1:5, function(ts) sd(impact_8[[ts]]))
impact_df <- data.frame(mean = c(mean_bEH_only, mean_LA_only, mean_both),
                        sd = c(sd_bEH_only, sd_LA_only, sd_both),
                        ts = rep(c("B", "Bd", "H", "E", "A"), times = 3),
                        int_type = rep(c("bEH", "LA", "both"), each = 5))

#impact = list(impact_6, impact_7, impact_8)
#impact_df <- adply(1:3, 1, function(i) {
#    adply(1:5, 1,function(ts) {
#        data.frame(impact = impacts[[i]][[ts]],
#                   ts = rep(c("B", "Bd", "H", "E", "A")[ts], times = length(impacts[[i]])),
#                   int_type = rep(c("bEH", "LA", "both")[i],  each = length(impacts[[i]])))
#               })
#})

#str(impact_df)
#impact_df <- ddply(impact_df, .(ts, int_type), summarise, mean = mean(impact), sd = sd(impact))
ggplot(impact_df, aes(int_type, mean, col = ts)) + geom_point() + facet_wrap(ts~.) +
    geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd))
"""
