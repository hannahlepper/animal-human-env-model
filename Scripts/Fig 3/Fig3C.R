#Get env model df
TS <- c("Baseline", "Balanced", "Human\ndominated", "Animal\ndominated",  "Environment\ndominated")
df_unbounded <- lapply(1:5, function(ts) {
  df <- data.frame(bEH = as.numeric(bEH_unif[1:200000]), 
                   impact = fig3C_impacts_unbounded[[ts]][1:200000],
                   TS = TS[ts]) 
}) %>% bind_rows() %>%
  mutate(., pv_cut = cut(bEH, seq(0., 1., 0.01)),
    TS = factor(TS, levels = c("Baseline", "Balanced", "Human\ndominated", 
                                           "Environment\ndominated", "Animal\ndominated"), 
                            ordered = TRUE)) %>%
  ddply(., .(pv_cut, TS), summarise, 
    bEH = mean(bEH),
    n= length(impact),
    n_na = sum(is.na(impact)),
    n_nan = sum(is.nan(impact)),
    mean_impact = mean(impact),
    med_impact = median(impact), 
    sd_impact = sd(impact),
    #impact_sample = paste0(impact[1:3]),
    upper = quantile(impact, probs = .75), 
    lower = quantile(impact, probs = .25))

print(head(df_unbounded))
#Get original data
df_orig_unbounded = data.frame(bEH = 0, TS = "No environment",
  mean_impact = mean(impact_orig), 
  med_impact = median(impact_orig),
  lower = quantile(impact_orig, probs = .25), 
  upper = quantile(impact_orig, probs = .75))

#dummy data
#df_unbounded = data.frame(impact = rnorm(10000), TS = sample(letters[1:4], 10000, replace = TRUE), bEH = rbeta(10000, 1,1))


png("plots/Fig3c_smooth.png", width = 10, height = 8, units = "cm", res = 300)
#p <- plot(df_unbounded$bEH, df_unbounded$impact)
p <- ggplot(df_unbounded, aes(bEH, mean_impact)) + 
    #geom_pointrange(data = df_orig_unbounded, 
    #    mapping = aes(x = bEH, y = mean_impact, ymin = lower, ymax = upper), 
    #    colour = "black", size = .2) +
    geom_line(aes(col = TS)) +
    #stat_smooth(method = "loess") +
    viridis::scale_colour_viridis() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = TS), 
      show.legend = FALSE, alpha = 0.2) +
    theme_bw(base_size = 7) +
    labs(y = expression(paste(omega)),
         x = expression(paste(beta[EH])),        
         col = "Transmission\nscenario") 
    #theme(legend.background = element_blank())
print(p)
dev.off()