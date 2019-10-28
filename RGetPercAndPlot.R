library(plyr)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

get_perc_target <- function(lower_bin, p, dat) {
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

get_mean_var <- function(lower_bin, p, dat) {
  df <- data.frame(p, dat)
  names(df) <- c("x", "y", "dat")
  df <- mutate(df, x_cat = cut(df$x, c(lower_bin, max(df$x)), 
                               labels = as.character(lower_bin), 
                               orderd_result = TRUE)) %>%
    mutate(., y_cat = cut(df$y, c(lower_bin, max(df$y)),
                          labels = as.character(lower_bin),
                          ordered_result = TRUE)) %>%
    ddply(., .(x_cat, y_cat), summarise, mean = mean(dat), var = var(dat))
  return(df)
}

myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))

plot_heatmap <- function(df, axisnames, plottitle, filltitle) {
  ggplot(df, aes(x_cat, y_cat, fill = mean)) + 
    geom_tile() +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour= "black"),
          panel.grid = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          strip.background = element_rect(fill = NA, colour = NA),
          axis.text = element_text(size = 8)) +
    scale_fill_gradientn(colors = myPalette(100), limits = c(0,1)) +
    labs(x = axisnames[1], y = axisnames[2], title = plottitle) 
}

plot_heatmap_var <- function(df, axisnames, plottitle, filltitle) {
  ggplot(df, aes(x_cat, y_cat, fill = var)) + 
    geom_tile() +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour= "black"),
          panel.grid = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          strip.background = element_rect(fill = NA, colour = NA),
          axis.text = element_text(size = 8)) +
    scale_fill_gradientn(colors = myPalette(100), limits = c(0,1)) +
    labs(x = axisnames[1], y = axisnames[2], title = plottitle) 
}

#testing
p <- matrix(rlnorm(1000, meanlog = log(0.01)+2, sdlog = 2), ncol = 2)
dat <- rbinom(5000, 1, 0.5)
lower_bin <- seq(0., 1., 0.05)

df <- get_perc_target(lower_bin, p, dat)
p1 <- plot_heatmap(df, c("a", "b"), "c", "d")
svg("M:/test.svg", height = 12, width = 20)
grid.arrange(p1, p1, p1, p1, p1, p1, p1, p1, p1, p1, p1, p1, nrow = 4, ncol = 3)
dev.off()
