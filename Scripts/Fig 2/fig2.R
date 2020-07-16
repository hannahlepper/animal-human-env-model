#Pv 1 and 4.
#Impacts 3 and 4
dfs <- adply(1:5, 1, function(ts) {
    data.frame(param_val = c(Pv[[1]][[ts]][1:10000,11], Pv[[4]][[ts]][1:10000,2]),
               param = rep(c("bEH", "LA"), each = 10000),
               impact = c(impact_3[[ts]][1:10000], impact_4[[ts]][1:10000]))
})
head(dfs)
dfs = mutate(dfs, ts = rep(c("B", "Bd", "H", "E", "A"), each = 20000)) %>%
      mutate(., pv_cut = cut(dfs$param_val, c(0, 0.095, 0.105, 0.495, 0.505, 1.0)))
df_low <- dfs[dfs$pv_cut == "(0.095,0.105]",]
df_high <- dfs[dfs$pv_cut == "(0.495,0.505]",]

png("Plots/Fig2.png")
wrap_plots(ggplot(df_low, aes(ts, impact, col = param)) + geom_boxplot(),
      ggplot(df_high, aes(ts, impact, col = param)) + geom_boxplot())
dev.off()
