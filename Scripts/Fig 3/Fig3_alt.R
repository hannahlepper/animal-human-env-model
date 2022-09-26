TS <- c("Balanced", "Human\ndominated", "Animal\ndominated", "Environment\ndominated")

df_unbounded <- lapply(1:length(f3_impacts_unbounded), function(ts) {
    data.frame(TS = TS[ts], impact = f3_impacts_unbounded[[ts]], mod = "Unbounded")
}) %>%
    bind_rows()

df_bounded <- lapply(1:length(f3_impacts_bounded), function(ts) {
    data.frame(TS = TS[ts], impact = f3_impacts_bounded[[ts]], mod = "Bounded")
}) %>%
    bind_rows()

df_orig <- lapply(1:length(f3_impacts_original), function(ts) {
    data.frame(TS = TS[ts], impact = f3_impacts_original[[ts]], mod = "Original")
}) %>%
    bind_rows()


df <- bind_rows(df_unbounded, df_bounded, df_orig) %>%
    ddply(., .(TS, mod), summarise, 
        mean_impact = mean(impact), 
        med_impact = median(impact),
        lower = quantile(impact, .25), 
        upper = quantile(impact, .75), 
        sd = sd(impact)) %>%
    mutate(., TS = factor(TS, levels = c("Balanced", "Human\ndominated", 
                                           "Environment\ndominated", "Animal\ndominated"), 
                            ordered = TRUE))

print(df)


png("plots/Fig3c_alt.png", width = 15, height = 8, units = "cm", res = 300)
p <- ggplot(df, aes(TS, mean_impact, col = mod)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), , shape = 1,
        position = position_dodge(width = 0.2)) +
    labs(y = expression("Intervention impact (" ~ omega ~ ")"), x = "", col = "Environment type") +
    viridis::scale_colour_viridis(discrete = TRUE) +
    theme_bw() +
    theme(legend.position = c(0.12, 0.8),
            legend.background = element_blank()) 
print(p)
dev.off()