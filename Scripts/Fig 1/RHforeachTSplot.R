source('RGetPercAndPlot.R')
mean_impact <- lapply(c(15,13,1), function(p) {
        lapply(1:5, function(x) {
            get_perc_target(lower_bin, Pv[[1]][[x]][, c(11,p)], n_target[[x]])
        })
})
names(mean_impact) = c("muE", "muH", "LH")

#plotting!
bEH_mean_impact_plots <- lapply(c("muE", "muH", "LH"), function(p) {
    lapply(1:5, function(x) {
        plot_heatmap(mean_impact[[p]][[x]], c('bEH', p),
                     paste(c("% Target Achieved, ts = ", c("B", "Bd", "H", "E", "A")[x]), collapse = ""),
        limits = c(0, 0.75))
        })
    })

#svg('M:/Project folders/Model env compartment/Plots/ptaplot.svg', height=20, width = 20)
wrap_plots(unlist(bEH_mean_impact_plots, recursive = FALSE), ncol = 5)
#dev.off(