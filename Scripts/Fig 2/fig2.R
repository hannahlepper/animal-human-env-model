library(data.table)
library(plotmath)

#Pv 1 and 4.
#Impacts 3 and 4
TS <- c("Baseline", "Balanced", "Human\ndominated", 
          "Environment\ndominated", "Animal\ndominated")
dfs <- adply(1:5, 1, function(ts) {
    data.frame(param_val = c(P[[1]][[ts]][,11], P[[2]][[ts]][,2]),
               param = rep(c("beta[EH]", "Lambda[A]"), each = dim(impact_3[[ts]])[1]),
               impact = c(impact_3[[ts]], impact_4[[ts]]),
               ts = rep(TS[ts], dim(impact_4[[ts]])[1] * 2))
}) %>%
      data.table()

df_low <- dfs[param_val >= 0.095 & param_val <= 0.105]
df_high <- dfs[param_val >= 0.495 & param_val <= 0.505]

label_parse <- function(breaks) parse(text = breaks)

png("plots/boxplotFig2A", width = 15, height = 10, units = "cm", res = 300)
ggplot(df_low, aes(ts, impact, col = param)) + 
      geom_boxplot() +
      labs(x = "", y = "Impact") +
      scale_colour_discrete("Parameter targeted\nfor intervention",labels = label_parse) +
      theme_bw()+
      theme(legend.position = c(0.15, 0.85),
            legend.background = element_blank())
dev.off()

png("plots/boxplotFig2B", width = 15, height = 10, units = "cm", res = 300)
ggplot(df_high, aes(ts, impact, col = param)) + 
      geom_boxplot() +
      labs(x = "", y = "Impact") +
      scale_colour_discrete("Parameter targeted\nfor intervention",labels = label_parse) +
      theme_bw()+
      theme(legend.position = c(0.15, 0.85),
            legend.background = element_blank())
dev.off()
