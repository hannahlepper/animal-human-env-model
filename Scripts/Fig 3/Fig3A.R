#For testing - dummy data
Pv <- lapply(1:100, function(x) data.frame(matrix(runif(55, 0, 1), ncol = 11)))
impact_2 <- lapply(1:5, function(x) runif(5, 0, 0.15))
impact_5 <- lapply(1:5, function(x) runif(5, 0, 0.15))

#myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
myP <- colorRampPalette(c("#070066","#1100FA","#BAD3F7", "#F6FCBD", "#FFC803"))

plot_heatmap2 <- function(df, title) {
  ggplot(df, aes(x, y, z = impact)) +
    stat_summary_2d(binwidth = 0.1) +
    scale_fill_gradientn(colors = myP(100), limits = c(0, 0.15)) +
    labs(x = expression(paste(beta[EH])), 
         y = expression(paste(Lambda[A], " (Pre-intervention)")),
         fill = expression(paste(omega[A])),
         title = title) +
    theme_bw()
}

#Data for heatmaps needed
#impacts_highbHA <- d_impact
#impacts_lowbHA <- d_impact
impacts_highbHA <- impact_2
impacts_lowbHA <- impact_5

TS <- c("Baseline", "Balanced", "Human-dominated", "Environment-dominated", "Animal-dominated")
df_list_highbHA <- lapply(1:5, function(ts) {
 paramdf <- data.frame(x = as.numeric(Pv[[10+ts]][,11]), y = as.numeric(Pv[[10+ts]][,2])) 
  df <- data.frame(paramdf, impacts_highbHA[ts]) %>%
    setNames(., c("x", "y", "impact"))  
})
df_list_lowbHA <- lapply(1:5, function(ts) {
  paramdf <- data.frame(x = as.numeric(Pv[[40+ts]][,11]), y = as.numeric(Pv[[40+ts]][,2])) 
  df <- data.frame(paramdf, impacts_lowbHA[ts]) %>%
    setNames(., c("x", "y", "impact")) 
})

#Get legend
legenddf <- df_list_highbHA[[1]]
legend_plot <- plot_heatmap2(legenddf, paste0(TS[1], " transmission scenario"))
legend_to_plot <- cowplot::get_legend(legend_plot)

plot_list_highbHA <- lapply(1:6, function(ts) {
  if (ts < 6) {
    df <- df_list_highbHA[[ts]]
    gplot <- plot_heatmap2(df, paste0(TS[ts]))
    gplot <- gplot + theme(legend.position = "none")
    return(gplot)
  } else {
    return(legend_to_plot)
  }
})

plot_list_lowbHA <- lapply(1:6, function(ts) {
  if (ts < 6) {
    df <- df_list_lowbHA[[ts]]
    gplot <- plot_heatmap2(df, paste0(TS[ts]))
    gplot <- gplot + theme(legend.position = "none")
    return(gplot)
  } else {
    return(legend_to_plot)
  }
})

png("plots/AppendixFig2highbHA.png", width = 18, height = 23, units = "cm", res = 300)
do.call(grid.arrange, plot_list_highbHA)
dev.off()

png("plots/AppendixFig2lowbHA.png", width = 18, height = 23, units = "cm", res = 300)
do.call(grid.arrange, plot_list_lowbHA)
dev.off()

png("plots/Fig3.png", width = 19, height = 8, units = "cm", res = 300)
cowplot::plot_grid(plot_list_highbHA[[2]], plot_list_lowbHA[[2]], legend_to_plot, 
  ncol = 3, rel_widths = c(0.4, 0.4, 0.2))
dev.off()

lapply(1:5, function(ts) {

  fn <- paste(c("plots/heatmap_highbHA_", TS[ts], ".png"), collapse = "")
  png(fn, width = 17, height = 13, units = "cm", res = 300)
  print(plot_list_highbHA[[ts]])
  dev.off()

  fn <- paste(c("plots/heatmap_lowbHA_", TS[ts], ".png"), collapse = "")
  png(fn, width = 17, height = 13, units = "cm", res = 300)
  print(plot_list_lowbHA[[ts]])
  dev.off()

})

