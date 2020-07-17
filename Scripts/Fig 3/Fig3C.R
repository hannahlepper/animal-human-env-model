impact_1 = d_impact
dfs <- data.frame(bEH_val = unlist(lapply(1:5, function(elem) Pv[[1]][[elem]][[11]])),
                  impact = unlist(lapply(1:5, function(elem) impact_1[[elem]])))

str(dfs)

impact_9 = d_impact[[1]]
dfs_orig <- data.frame(bEH_val = rep(0, length(impact_1[[1]])),
                       impact = impact_9)
dfs2 = rbind(dfs, dfs_orig)
dfs2 = mutate(dfs2, ts = rep(c('Baseline', 'Balanced', 'Human\ndominated', 'Environment\ndominated', 'Animal\ndominated', 'Original: no\nenvironment'), each = length(impact_1[[1]]))) %>%
      mutate(., pv_cut = cut(dfs2$bEH_val, seq(0., 1., 0.01)))

prop_sd <- function(p, N) {
  if(p > 0 & N > 0) {
    sqrt((p*(1-p)))/sqrt(N)
  } else {
    0
  }
}

dfs_summary <- ddply(dfs2, .(ts, pv_cut), summarise, 
                     mean_impact = mean(impact), 
                     sd_impact = prop_sd(mean(impact), length(impact)), 
                     bEH = mean(bEH_val))


png("plots/Fig3c_dotplot.png", width = 15, height = 8, units = "cm", res = 300)
p <- ggplot(dfs_summary, aes(bEH, mean_impact, col = ts)) + geom_point(position = 'jitter') +
    geom_errorbar(aes(ymin = mean_impact - sd_impact, ymax = mean_impact + sd_impact)) +
    geom_line() +
    theme_bw(base_size = 7) +
    labs(y = parse(text = expression(omega[A])), 
         x = parse(text = expression(beta[EH])), 
         col = "Transmission\nscenario")+
    theme(legend.position = c(0.9, 0.65),
          legend.background = element_blank())
print(p)
dev.off()