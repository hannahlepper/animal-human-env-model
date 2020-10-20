

#myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
myP <- colorRampPalette(c("#070066","#1100FA","#428af5","#BAD3F7", "#F6FCBD", "#FFC803"),
                        bias = 1.)

plot_heatmap2 <- function(df, title) {
  ggplot(df, aes(x, y, z = impact)) +
    stat_summary_2d() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_gradientn(colors = myP(100), limits = c(0, .127)) +
    labs(title = title) 
}

#Data for heatmaps needed
#impacts_highbHA <- d_impact
#impacts_lowbHA <- d_impact
impacts_highbHA <- impact_2
impacts_lowbHA <- impact_5

TS <- c("Baseline", "Balanced", "Humandominated", "Environmentdominated", "Animaldominated")
plot_list_highbHA <- lapply(1:5, function(ts) {
  #mean_impact = get_mean(lower_bin, Pv[[15+ts]][,c(11,2)], impacts_highbHA[ts])
  df <- data.frame(Pv[[10+ts]][,c(11,2)], impact = impacts_highbHA[ts]) %>%
    setNames(., c("x", "y", "impact"))
  #plot_heatmap(mean_impact, paste(c(TS[ts], " high bHA"), collapse = ""))
  plot_heatmap2(df, paste(c(TS[ts], " high bHA"), collapse = ""))
})

plot_list_lowbHA <- lapply(1:5, function(ts) {
  #mean_impact = get_mean(lower_bin, Pv[[45+ts]][,c(11,2)], impacts_lowbHA[ts])
  df <- data.frame(Pv[[40+ts]][,c(11,2)], impact = impacts_lowbHA[ts]) %>%
    setNames(., c("x", "y", "impact"))
  #plot_heatmap(mean_impact, paste(c(TS[ts], " low bHA"), collapse = ""))
  plot_heatmap2(df, paste(c(TS[ts], " low bHA"), collapse = ""))
})

lapply(1:5, function(ts) {

  fn <- paste(c("plots/heatmap_highbHA_", TS[ts], ".png"), collapse = "")
  png(fn, width = 8, height = 7, units = "cm", res = 300)
  print(plot_list_highbHA[[ts]])
  dev.off()

  fn <- paste(c("plots/heatmap_lowbHA_", TS[ts], ".png"), collapse = "")
  png(fn, width = 8, height = 7, units = "cm", res = 300)
  print(plot_list_lowbHA[[ts]])
  dev.off()

})