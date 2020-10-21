

#myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
myP <- colorRampPalette(c("#070066","#1100FA","#428af5","#BAD3F7", "#F6FCBD", "#FFC803"),
                        bias = 1.)

plot_heatmap2 <- function(df, title) {
  ggplot(df, aes(x, y, z = impact)) +
    stat_summary_2d(binwidth = 0.1) +
    scale_fill_gradientn(colors = myP(100), limits = c(0, .127)) +
    labs(x = expression(paste("Transmission from environment to humans (", beta[EH], ")")), 
         y = expression(paste("Pre-intervention antibiotic consumption in animals (", Lambda[A], ")")),
         title = title) +
    theme_bw()
}

#Data for heatmaps needed
#impacts_highbHA <- d_impact
#impacts_lowbHA <- d_impact
impacts_highbHA <- impact_2
impacts_lowbHA <- impact_5

TS <- c("Baseline", "Balanced", "Humandominated", "Environmentdominated", "Animaldominated")
df_list_highbHA <- lapply(1:5, function(ts) {
 paramdf <- data.frame(x = as.numeric(Pv[[40+ts]][,11]), y = as.numeric(Pv[[40+ts]][,2])) 
  df <- data.frame(paramdf, impacts_highbHA[ts]) %>%
    setNames(., c("x", "y", "impact"))  
})
df_list_lowbHA <- lapply(1:5, function(ts) {
  paramdf <- data.frame(x = as.numeric(Pv[[40+ts]][,11]), y = as.numeric(Pv[[40+ts]][,2])) 
  df <- data.frame(paramdf, impacts_lowbHA[ts]) %>%
    setNames(., c("x", "y", "impact")) 
})

plot_list_highbHA <- lapply(df_list_highbHA, function(df) {
  plot_heatmap2(df, paste(c(TS[ts], " high bHA"), collapse = ""))
})

plot_list_lowbHA <- lapply(df_list_lowbHA, function(df) {
  plot_heatmap2(df, paste(c(TS[ts], " low bHA"), collapse = ""))
})


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