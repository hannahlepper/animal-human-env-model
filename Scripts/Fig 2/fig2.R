library(data.table)

#dummy data for testing
#impact_3 = d_impact
#impact_4 = d_impact
#Pv = P

#Pv 1 and 4.
#Impacts 3 and 4
TS <- c("Balanced", "Human\ndominated", 
          "Animal\ndominated", "Environment\ndominated")
dfs1 <- adply(1:4, 1, function(ts) {
    data.frame(param = rep(c("beta[EH]", "Lambda[A]"), 
                           each = length(bEH_impact_unbounded[[ts]])),
               impact = c(bEH_impact_unbounded[[ts]], LA_impact_unbounded[[ts]]),
               ts = rep(TS[ts], length(LA_impact_unbounded[[ts]]) * 2))
}) %>%
      data.table() %>%
      mutate(., ts = factor(ts, levels = c("Balanced", "Human\ndominated", 
                                           "Environment\ndominated", "Animal\ndominated"), 
                            ordered = TRUE),
                model = "Unbounded environment")

dfs2 <- adply(1:4, 1, function(ts) {
    data.frame(param = rep(c("beta[EH]", "Lambda[A]"), 
                           each = length(bEH_impact_bounded[[ts]])),
               impact = c(bEH_impact_bounded[[ts]], LA_impact_bounded[[ts]]),
               ts = rep(TS[ts], length(LA_impact_bounded[[ts]]) * 2))
}) %>%
      data.table() %>%
      mutate(., ts = factor(ts, levels = c("Balanced", "Human\ndominated", 
                                           "Environment\ndominated", "Animal\ndominated"), 
                            ordered = TRUE),
                 model = "Bounded environment")

dfs <- bind_rows(dfs1, dfs2)
label_parse <- function(breaks) parse(text = breaks)

print(dim(dfs))
print(unique(dfs$model))
print(head(dfs2))

dfs_summary <- ddply(dfs, .(ts, param, model), summarise, 
      mean_impact = mean(impact), lwr = quantile(impact, 0.25), 
      upr = quantile(impact, 0.75))

png("plots/boxplotFig2A.png", width = 20, height = 12, units = "cm", res = 300)
p <- ggplot(dfs_summary, aes(ts, mean_impact, col = param)) + 
      geom_pointrange(aes(ymin = lwr, ymax = upr), shape = 1,
                     position = position_dodge(width = 0.2)) +
      labs(x = "", y = expression(omega)) +
      scale_colour_discrete("Parameter targeted\nfor intervention", labels = label_parse) +
      theme_bw()+
      theme(legend.position = c(0.1, 0.85),
            legend.background = element_blank()) +
      facet_wrap(~model, ncol = 2)
print(p)
dev.off()

