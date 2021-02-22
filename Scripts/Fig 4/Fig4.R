myP <- colorRampPalette(c("#070066","#1100FA","#BAD3F7", "#F6FCBD", "#FFC803"))
myP <- colorRampPalette(c("#1100FA","#FFFFFF", "#FFC803"))


plot_heatmap2 <- function(df, title) {
  ggplot(df, aes(x, y, z = impact)) +
    stat_summary_2d(binwidth = 0.1) +
    scale_fill_gradientn(colors = myP(100), limits = c(0, 0.5)) +
    labs(x = expression(paste("Transmission from environment to humans (", beta[EH], ")")), 
         y = expression(paste("Pre-intervention antibiotic consumption in animals (", Lambda[A], ")")),
         title = title) +
    theme_bw()
}

df_list <- lapply(1:5, function(ts) {
  paramdf <- data.frame(x = as.numeric(Pv[[10+ts]][,11]), y = as.numeric(Pv[[10+ts]][,2])) 
  df <- data.frame(paramdf, impact_10[ts]) %>%
    setNames(., c("x", "y", "impact")) 
})

plot_list <- lapply(1:5, function(ts) {
  df <- df_list[[ts]]
  plot_heatmap2(df, paste(c(TS[ts], " double intervention"), collapse = ""))
})

lapply(1:5, function(ts) {

  fn <- paste(c("plots/heatmap_doubleintervention_", TS[ts], ".png"), collapse = "")
  png(fn, width = 17, height = 13, units = "cm", res = 300)
  print(plot_list[[ts]])
  dev.off()

})