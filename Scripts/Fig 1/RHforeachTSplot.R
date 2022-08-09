#sort out julia's list
R <- unlist(R)
R2 <- unlist(R2)
#For testing - dummy data
#R = rbeta(3 * 5, 1, 1)
Rdf <- data.frame(
        R = R, model = "Unbounded\nenvironment",
        R_type = factor(rep(c("RH", "RA", "RE"), times = 4)),
        TS = rep(c("Balanced", "Environment\ndominated", "Animal\ndominated", "Human\ndominated"), each = 3))

Rdf2 <- data.frame(
        R = R2, model = "Bounded\nenvironment",
        R_type = factor(rep(c("RH", "RA", "RE"), times = 4)),
        TS = rep(c("Balanced", "Environment\ndominated", "Animal\ndominated", "Human\ndominated"), each = 3))

Rdf <- bind_rows(Rdf, Rdf2) %>%
    subset(., (R_type == "RE")) %>%
    mutate(., TS = factor(TS, 
                        levels = c("Balanced", "Human\ndominated", 
                                    "Environment\ndominated", "Animal\ndominated"),
                        ordered = TRUE),
             R_type = factor(R_type, levels = c("RH", "RA", "RE"),
                        labels = c("RH", "Proportion of animals carrying\n resistant bacteria",
                                    "Proportion of environmental units\n carrying resistant bacteria")))

png("plots/barplotFig1B.png", width = 9, height = 8, units = "cm", res = 300)
p <- ggplot(Rdf, aes(TS, R, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    #scale_fill_discrete("", labels = parse(text = levels(Rdf$R_type))) +
    labs(x = "", y = "Proportion of environmental units\n carrying resistant bacteria", fill = "Model") +
    #viridis::scale_fill_viridis() +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
            legend.position = c(0.15, 0.85),
            legend.background = element_blank())
    #facet_wrap(~R_type, ncol = 1)
print(p)
dev.off()
