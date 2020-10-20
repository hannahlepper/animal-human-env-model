#impact_1 = d_impact
dfs <- data.frame(bEH_val = unlist(lapply(2:5, function(elem) Pv[[elem]][,11])),
                  impact = unlist(lapply(2:5, function(elem) impact_1[[elem]])))

str(dfs)

#impact_9 = d_impact[[1]]
dfs_orig <- data.frame(bEH = rep(0, length(impact_1[[1]])),
                       impact = impact_9) %>%
               summarise(., med_impact = median(impact),
                         lower = quantile(impact, .25),
                         upper = quantile(impact, .75), 
                         bEH = 0)
#dfs2 = rbind(dfs, dfs_orig)
dfs2 = dfs
dfs2 = mutate(dfs2, ts = rep(c('Balanced', 'Human\ndominated', 'Environment\ndominated', 'Animal\ndominated'), each = length(impact_1[[1]]))) %>%
       mutate(., pv_cut = cut(dfs2$bEH_val, seq(0., 1., 0.01)))


dfs_summary <- ddply(dfs2, .(ts, pv_cut), summarise, 
                     med_impact = median(impact), 
                     lower = quantile(impact, .25),
                     upper = quantile(impact, .75), 
                     bEH = mean(bEH_val))


png("plots/Fig3c_dotplot.png", width = 15, height = 8, units = "cm", res = 300)
p <- ggplot(dfs_summary, aes(bEH, med_impact, col = ts)) + 
    geom_pointrange(data = dfs_orig, mapping = aes(ymin = lower, ymax = upper), colour = "black", size = .2) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=ts), alpha=0.1, colour = NA) +
    theme_bw(base_size = 7) +
    #labs(y = parse(text = expression(omega[A])), 
    #     x = parse(text = expression(beta[EH])), 
     labs(y = "Impact of reducing antibiotic consumption in animals",
          x = "Transmission from environment to humans",        
          col = "Transmission\nscenario") +
    theme(legend.background = element_blank())
print(p)
dev.off()