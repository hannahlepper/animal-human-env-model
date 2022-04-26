#For testing - dummy data
#Pv <- lapply(1:100, function(x) data.frame(matrix(runif(55, 0, 1), ncol = 11)))
#impact_2 <- lapply(1:5, function(x) runif(5, 0, 0.15))
#impact_5 <- lapply(1:5, function(x) runif(5, 0, 0.15))

#myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
#myP <- colorRampPalette(c("#070066","#1100FA","#BAD3F7", "#F6FCBD", "#FFC803"))
library(viridis)
plot_heatmap2 <- function(df, title, lim, breaks) {
  ggplot(df, aes(x, y, z = impact)) +
    stat_summary_2d(binwidth = 0.1) +
    #scale_fill_gradientn(colors = myP(100), limits = c(0, lim),
    #    breaks = breaks) +
    scale_fill_viridis(limits = c(0, lim)) +
    labs(x = expression(paste(beta[EH])), 
         y = expression(paste(Lambda[A], " (Pre-intervention)")),
         fill = expression(paste(omega[A])),
         title = title) +
    theme_bw()
}

#Data for appendix heatmaps needed

TS <- c("Baseline", "Balanced", "Human-dominated", "Animal-dominated",  "Environment-dominated")
df_list_unbounded <- lapply(1:5, function(ts) {
  paramdf <- data.frame(x = as.numeric(bEH_unif), y = as.numeric(LA_unif)) 
  df <- data.frame(paramdf, heatmap_impacts_unbounded[[ts]]) %>%
    setNames(., c("x", "y", "impact"))  
})
df_list_bounded <- lapply(1:5, function(ts) {
  paramdf <- data.frame(x = as.numeric(bEH_unif), y = as.numeric(LA_unif)) 
  df <- data.frame(paramdf, heatmap_impacts_bounded[[ts]]) %>%
    setNames(., c("x", "y", "impact")) 
})

#Get legend
legenddf <- df_list_unbounded[[1]]
legend_plot <- plot_heatmap2(legenddf, paste0(TS[1], " transmission scenario"), lim = 0.1, breaks = c(0, 0.05, 0.1, 0.1, 0.1))
legend_to_plot <- cowplot::get_legend(legend_plot)

plot_list_appendix_unbounded <- lapply(1:6, function(ts) {
  if (ts < 6) {
    df <- df_list_unbounded[[ts]]
    gplot <- plot_heatmap2(df, paste0(TS[ts]), lim = 0.1, breaks = c(0, 0.05, 0.1, 0.1, 0.1))
    gplot <- gplot + theme(legend.position = "none")
    return(gplot)
  } else {
    return(legend_to_plot)
  }
})

plot_list_appendix_bounded <- lapply(1:6, function(ts) {
  if (ts < 6) {
    df <- df_list_bounded[[ts]]
    gplot <- plot_heatmap2(df, paste0(TS[ts]), lim = 0.1, breaks = c(0, 0.05, 0.1, 0.1, 0.1))
    gplot <- gplot + theme(legend.position = "none")
    return(gplot)
  } else {
    return(legend_to_plot)
  }
})

png("plots/AppendixFig2Unbounded.png", width = 18, height = 23, units = "cm", res = 300)
do.call(grid.arrange, plot_list_appendix_bounded)
dev.off()

png("plots/AppendixFig2Bounded.png", width = 18, height = 23, units = "cm", res = 300)
do.call(grid.arrange, plot_list_appendix_unbounded)
dev.off()

#Fow figures for within the text
df_unbounded <- data.frame(x = as.numeric(bEH_unif), y = as.numeric(LA_unif),
                 impact = as.numeric(fig3A_impacts_unbounded)) %>%
    setNames(., c("x", "y", "impact"))
df_bounded <- data.frame(x = as.numeric(bEH_unif), y = as.numeric(LA_unif),
                 impact = as.numeric(fig3A_impacts_bounded)) %>%
    setNames(., c("x", "y", "impact"))

plot_lowbHA_unbounded <- plot_heatmap2(df_unbounded, "", .1, breaks = NA) + 
  theme(legend.position = "none") 
#plot_lowbHA_bounded <- plot_heatmap2(df_bounded, "", .1, breaks = NA) + 
#  theme(legend.position = "none")


legend_plot <- plot_heatmap2(df_unbounded, "", lim = 0.1, breaks = NA)
legend_to_plot <- cowplot::get_legend(legend_plot)

png("plots/Fig3.png", width = 8, height = 8, units = "cm", res = 300)
print(plot_lowbHA_unbounded)
#cowplot::plot_grid(plot_list_highbHA[[2]], plot_list_lowbHA[[2]], legend_to_plot, 
#  ncol = 3, rel_widths = c(0.4, 0.4, 0.2))
dev.off()

png("plots/Fig3Legend.png", width = 5, height = 5, units = "cm", res = 300)
plot(legend_to_plot)
dev.off()

df_unbounded_summmary <- df_unbounded %>%
  mutate(., beh_cut = cut(x, breaks = seq(0, 1, 0.1)), 
            LA_cut = cut(y, breaks = seq(0, 1, 0.1))) %>%
  ddply(., .(beh_cut, LA_cut), summarise,
    mean_impact = mean(impact))
  
print(df_unbounded_summmary)