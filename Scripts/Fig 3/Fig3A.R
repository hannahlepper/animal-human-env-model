get_mean <- function(lower_bin, p, dat) {
  df <- data.frame(p, dat)
  names(df) <- c("x", "y", "dat")
  df <- mutate(df, x_cat = cut(df$x, c(lower_bin, max(df$x)),
                               labels = as.character(lower_bin),
                               orderd_result = TRUE)) %>%
    mutate(., y_cat = cut(df$y, c(lower_bin, max(df$y)),
                          labels = as.character(lower_bin),
                          ordered_result = TRUE)) %>%
    ddply(., .(x_cat, y_cat), summarise, mean = mean(dat))
  return(df)
}

#myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
myP <- colorRampPalette(c("#070066","#1100FA","#428af5","#BAD3F7", "#F6FCBD", "#FFC803"),
                        bias = 1.)

plot_heatmap <- function(df, title) {
  ggplot(df, aes(x_cat, y_cat, fill = mean)) +
    geom_tile() +
    theme_bw(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradientn(expression(omega[A]), colors = myP(100), limits = c(0, .15)) +
    ylab(expression(Lambda[A] ~ "(pre-intervention)")) +
    xlab(expression(beta[EH])) + 
    labs(title = title)
}

p1 <- plot_heatmap(get_mean(lower_bin, P[[1]][[1]][,c(11,2)], d_impact[[1]]), "b")

#Data for heatmaps needed
impacts_highbHA <- d_impact
impacts_lowbHA <- d_impact
impacts_highbHA <- impact_2
impacts_lowbHA <- impact_5

TS <- c("Baseline", "Balanced", "Human\ndominated", "Environment\ndominated", "Animal\ndominated")
plot_list_highbHA <- lapply(1:5, function(ts) {
  mean_impact = get_mean(lower_bin, P[[3]][[ts]][,c(11,2)], impacts_highbHA[ts])
  plot_heatmap(mean_impact, paste(c(TS[ts], " high bHA"), collapse = ""))
})

plot_list_lowbHA <- lapply(1:5, function(ts) {
  mean_impact = get_mean(lower_bin, P[[3]][[ts]][,c(11,2)], impacts_lowbHA[ts])
  plot_heatmap(mean_impact, paste(c(TS[ts], " low bHA")))
})

TSfn <- c("Baseline", "Balanced", "Humandominated", "Environmentdominated", "Animaldominated")
lapply(1:5, function(ts) {

  fn <- paste(c("plots/heatmap_highbHA_", TSfn[ts]), collapse = "")
  png(fn, width = 8, height = 7, units = "cm", res = 300)
  print(plot_list_highbHA[[ts]])
  dev.off()

  fn <- paste(c("plots/heatmap_lowbHA_", TSfn[ts]), collapse = "")
  png(fn, width = 8, height = 7, units = "cm", res = 300)
  print(plot_list_lowbHA[[ts]])
  dev.off()

})