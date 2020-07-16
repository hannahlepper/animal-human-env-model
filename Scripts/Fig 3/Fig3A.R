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

plot_heatmap <- function(df) {
  ggplot(df, aes(x_cat, y_cat, fill = mean)) +
    geom_tile() +
    theme_bw(base_size = 5) +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradientn(expression(omega[A]), colors = myP(100), limits = c(0, 1)) +
    ylab(expression(Lambda[A] ~ "(pre-intervention")) +
    xlab(expression(beta[EH]))
}

p1 <- plot_heatmap(get_mean(lower_bin, P[[1]][[1]][,c(11,2)], d_impact[[1]]))
